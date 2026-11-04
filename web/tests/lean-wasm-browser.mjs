// Actual static page -> Lean WASM kernel -> C++ WASM, with no compiler API.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { createServer } from 'node:http';
import { mkdtemp, readFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { extname, resolve, sep } from 'node:path';
import { chromium } from 'playwright';
import { examples } from '../lean-study-examples.js';
const root = resolve(process.argv[2] || new URL('..', import.meta.url).pathname);
const prefix = '/sample-size-calc/';
const temporary = await mkdtemp(resolve(tmpdir(),'ssc-lean-browser-'));
const requests=[], errors=[], warnings=[];
const types={'.html':'text/html','.js':'text/javascript','.mjs':'text/javascript','.wasm':'application/wasm','.json':'application/json','.css':'text/css'};
const server=createServer(async(req,res)=>{
  const path=new URL(req.url,'http://localhost').pathname;
  requests.push({path,method:req.method});
  if(path.includes('/api/')) {res.writeHead(503).end('Backend disabled');return;}
  if(!path.startsWith(prefix)) {res.writeHead(404).end();return;}
  const file=resolve(root,path.slice(prefix.length)||'index.html');
  if(!file.startsWith(root+sep)) {res.writeHead(404).end();return;}
  try {
    const bytes=await readFile(file);
    // Match static Pages hosting: a subpath with no isolation headers.
    res.writeHead(200,{'Content-Type':types[extname(file)]||'application/octet-stream'}).end(bytes);
  } catch {res.writeHead(404).end();}
});
await new Promise(r=>server.listen(0,'127.0.0.1',r));
const browser=await chromium.launch({executablePath:process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH||undefined,args:['--disable-dev-shm-usage',...JSON.parse(process.env.SSC_CHROMIUM_ARGS||'[]')]});
const context=await browser.newContext({viewport:{width:1400,height:1000}});
const page=await context.newPage();
page.on('pageerror',e=>errors.push(e.message));
page.on('console',m=>{if(m.type()==='error') warnings.push(m.text());});
const base=`http://127.0.0.1:${server.address().port}${prefix}`;
let lastProgress='';
const progressTimer=setInterval(async()=>{
  try { const text=await page.locator('#lean-connection-status').textContent(); if(text!==lastProgress) {console.log('Compiler:',text);lastProgress=text;} } catch {}
},10000);
const compile=async source=>{
  if(source) await page.locator('#lean-source').fill(source);
  await page.locator('#lean-compile').click();
  await page.waitForFunction(()=>document.getElementById('lean-status').textContent==='Lean compilation complete.'||!document.getElementById('lean-error').hidden,null,{timeout:120000});
  assert.equal(await page.locator('#lean-error').isVisible(),false,await page.locator('#lean-error').textContent());
};
const ready=async()=>{
  await page.waitForFunction(()=>!document.getElementById('lean-compile').disabled || !document.getElementById('lean-setup').hidden,null,{timeout:300000});
  assert.equal(await page.locator('#lean-compile').isDisabled(),false,await page.locator('#lean-connection-status').textContent());
};
const download=async(id,name)=>{
  const waiting=page.waitForEvent('download');
  await page.locator(id).click();
  const result=await waiting;
  await result.saveAs(resolve(temporary,name));
  return resolve(temporary,name);
};
const reject=async(source,pattern)=>{
  await page.locator('#lean-source').fill(source);
  await page.locator('#lean-compile').click();
  await page.waitForFunction(()=>!document.getElementById('lean-error').hidden,null,{timeout:120000});
  assert.match(await page.locator('#lean-error').textContent(),pattern);
  console.log('Rejected invalid source:',pattern.toString());
};
try {
  await page.goto(base+'lean.html');
  await ready();
  assert.equal(await page.evaluate(()=>crossOriginIsolated),true,'service worker supplies browser isolation on a static host');
  assert.match(await page.locator('#lean-connection-status').textContent(),/in this browser/);
  console.log('Actual Lean WASM initialized.');
  await compile();
  assert.match(await page.locator('#lean-solution').textContent(),/51 participants per arm/);
  assert.match(await page.locator('#lean-theorem-type').textContent(),/variance study = 8 \/ 5/);
  assert.match(await page.locator('#lean-study-numbers').textContent(),/79.8176%/);
  const checked=JSON.parse(await readFile(await download('#lean-result-save','result.json'),'utf8'));
  await page.locator('#lean-cpp-details').evaluate(element=>{element.open=true;});
  const cpp=await download('#lean-cpp-save','study.cpp');
  const executable=resolve(temporary,'study');
  execFileSync('c++',['-std=c++17','-O2',cpp,'-o',executable],{timeout:60000});
  const native=JSON.parse(execFileSync(executable,{encoding:'utf8',timeout:10000}));
  const numerical=checked.results[0].evaluation;
  for(const [key,value] of Object.entries(numerical)) {
    if(typeof value==='number') assert.ok(Math.abs(native[key]-value)<=1e-12*Math.max(1,Math.abs(value)),key);
    else assert.deepEqual(native[key],value,key);
  }
  console.log('Downloaded C++ matches browser WASM.');
  // Independently cross-checked using SciPy normal tails.
  for(const [name,n] of Object.entries({'cluster-trial':18,'unequal-allocation':85,'targeted-ate':504,'custom-influence':252})) {
    await compile(examples[name]);
    assert.match(await page.locator('#lean-solution').textContent(),new RegExp(`n = ${n} `));
    console.log('Study passed:',name);
  }
  await compile(examples['repeated-measures']+'\nsolve power at100 using study at 100\nsolve effect effect100 using study at 100');
  assert.equal(await page.locator('#lean-query option').count(),3);
  await page.locator('#lean-query').selectOption('1'); assert.match(await page.locator('#lean-solution').textContent(),/Power =/);
  await page.locator('#lean-query').selectOption('2'); assert.match(await page.locator('#lean-solution').textContent(),/Detectable/);
  // The source is elaborated. A new function and changed effect cannot reuse a fixture.
  const changed='def plannedEffect (x : ℚ) : ℚ := x / 2\n'+examples['repeated-measures'].replace('means := [0, 0, 0, 1/2]','means := [0, 0, 0, plannedEffect 1/2]');
  assert.notEqual(changed,examples['repeated-measures']);
  await compile(changed);
  assert.match(await page.locator('#lean-solution').textContent(),/201 participants per arm/);
  await context.setOffline(true);
  // A fresh worker must load the complete compiler from persistent caches.
  await page.reload();
  await ready();
  await compile(examples['repeated-measures']);
  assert.match(await page.locator('#lean-solution').textContent(),/51 participants per arm/);
  console.log('New Lean definitions and offline reload/compilation passed.');
  await context.setOffline(false);
  await reject('def invalid : Nat := "wrong type"\n'+examples['repeated-measures'],/type mismatch/i);
  await reject(examples['cluster-trial'].replace('def icc : ℚ := 1/20','def icc : ℚ := 2'),/kernel|false|unapproved axiom/i);
  await reject(examples['targeted-ate'].replaceAll('propensity := 1/2','propensity := 0'),/positivity/i);
  await reject('theorem invalid : False := by sorry\n'+examples['repeated-measures'],/sorry/i);
  console.log('Checking compiler cancellation.');
  await page.locator('#lean-source').fill('#eval IO.sleep 60000\n'+examples['repeated-measures']);
  await page.locator('#lean-compile').click();
  await page.waitForTimeout(500);
  await page.locator('#lean-stop').click();
  await page.waitForFunction(()=>document.getElementById('lean-status').textContent==='Stopped.');
  await ready();
  await compile(examples['repeated-measures']);
  assert.deepEqual(requests.filter(r=>r.path.includes('/api/')),[]);
  assert.deepEqual(errors,[]);
  if(process.env.SCREENSHOTS) await page.screenshot({path:resolve(process.env.SCREENSHOTS,'lean-wasm.png'),fullPage:true});
  console.log('Browser Lean WASM: five statistical designs, altered source, all query targets, errors, cancellation and offline compilation passed. No compiler API requests.');
} catch(error) {
  console.error('Browser status:',await page.locator('#lean-connection-status').textContent());
  console.error('Browser errors:',errors,warnings.slice(-8));
  throw error;
} finally {clearInterval(progressTimer);await context.close();await browser.close();await new Promise(r=>server.close(r));await rm(temporary,{recursive:true,force:true});}
