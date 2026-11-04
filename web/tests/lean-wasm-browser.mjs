// Actual static page -> Lean WASM kernel and interpreter, with no compiler API.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { mkdtemp, readFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { extname, resolve, sep } from 'node:path';
import { chromium } from 'playwright';
import { examples } from '../lean-proof-examples.js';
import { checkCalculations } from './lean-calculations.mjs';
import { checkDerivations } from './lean-derivations.mjs';
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
const report = async name => JSON.parse(await readFile(await download('#lean-result-save', name),'utf8'));
try {
  await page.goto(base+'lean.html');
  await ready();
  assert.equal(await page.evaluate(()=>crossOriginIsolated),true,'service worker supplies browser isolation on a static host');
  assert.match(await page.locator('#lean-connection-status').textContent(),/in this browser/);
  console.log('Actual Lean WASM initialized.');
  await compile(examples['exact-trial']);
  assert.match(await page.locator('#lean-evaluation').textContent(),/some 25/);
  let checked = await report('trial.json');
  assert.equal(checked.compiler, 'lean-wasm-proofs-v2');
  const minimum = checked.results.find(d => d.name === 'Trial.minimum_sample_size');
  assert.ok(minimum);
  assert.match(minimum.statement, /adequate 25/);
  assert.ok(minimum.steps.some(s => s.before.includes('certificate') && s.after === 'no goals'));
  assert.ok(minimum.dependencies.includes('SSC.SearchResult.sound'));
  assert.ok(!minimum.axioms.some(a => /sorry|ofReduce|trustCompiler/.test(a)));
  const original = await download('#lean-proof-save','checked.lean');
  assert.equal(await readFile(original,'utf8'), examples['exact-trial']);
  console.log('Exact trial: witness, feasibility and minimality checked; proof trace exported.');

  for (const name of ['variance-propagation', 'certified-bounds']) {
    await compile(examples[name]);
    checked = await report(name+'.json');
    assert.ok(checked.results.some(d => d.kind === 'theorem' && d.parameters.length > 2 && d.steps.length > 0));
    if (name === 'certified-bounds') {
      const theorem = checked.results.find(d => d.name === 'Bounds.minimum_from_bounds');
      assert.ok(theorem.parameters.some(p => p.name === 'hEarlier' && p.hypothesis));
      assert.ok(theorem.steps.some(s => s.tactic === 'grind' && s.after === 'no goals'));
    } else assert.match(await page.locator('#lean-evaluation').textContent(),/9.*16/);
    console.log('Ordinary Lean module passed:',name);
  }

  // New, non-monotone predicates, structures and namespaces use the same compiler.
  const arbitrary = `import SSC
namespace UserModel
structure Design where
  permitted : Nat → Bool
def design : Design := ⟨fun n => n == 7 || n == 12⟩
def admissible (n : Nat) : Prop := design.permitted n = true
instance (n : Nat) : Decidable (admissible n) := inferInstanceAs (Decidable (_ = _))
def answer := SSC.leastBelow admissible 20
def exhausted := SSC.leastBelow admissible 7
#eval answer.value?
#eval exhausted.value?
theorem witness : answer.value? = some 7 := by decide +kernel
theorem empty_domain : exhausted.value? = none := by decide +kernel
theorem certificate : admissible 7 ∧ ∀ m, m < 7 → ¬ admissible m := by
  have h := answer.sound
  rw [witness] at h
  exact ⟨h.2.1, h.2.2⟩
end UserModel`;
  await compile(arbitrary);
  checked = await report('arbitrary.json');
  assert.ok(checked.results.some(d => d.name === 'UserModel.certificate'));
  assert.ok(!checked.results.some(d => d.name.startsWith('Trial.')), 'environments do not leak across runs');
  assert.match(await page.locator('#lean-evaluation').textContent(),/some 7[\s\S]*none/);
  console.log('New predicate: non-monotone minimum and exhausted domain certified.');

  // Alter the actual trial equation. Its old conclusion must no longer certify.
  const changed = examples['exact-trial'].replace('weights 3 10 n', 'weights 4 10 n');
  assert.notEqual(changed, examples['exact-trial']);
  await reject(changed, /false|kernel|Unproved dependency|unsolved/i);
  // Evaluate and certify the changed definition, without a prescribed answer.
  await compile(changed.slice(0, changed.indexOf('theorem selected_value')) + '\nend Trial');
  assert.doesNotMatch(await page.locator('#lean-evaluation').textContent(), /some 25/);
  console.log('Changed statistical equation invalidates the original proof.');

  await compile('import SSC\n#eval (List.range 8).sum\n#check Nat');
  assert.match(await page.locator('#lean-evaluation').textContent(),/28/);
  checked = await report('commands.json'); assert.equal(checked.results.length, 0);
  await reject('import SSC\ndef invalid : Nat := "wrong type"',/source.lean:2:.*type mismatch/i);
  await reject('import SSC\ntheorem invalid : False := by sorry',/sorry/i);
  await reject('import SSC\naxiom missing : False\ntheorem invalid : False := missing',/Unproved dependency.*missing/i);
  await reject('import SSC\ntheorem unproved (P : Nat → Prop) : P 9 := by grind',/unsolved|failed/i);
  await reject('import SSC\ntheorem native : 7 < 8 := by native_decide',/Unproved dependency|native_decide|ofReduce/i);
  console.log('False claims, missing hypotheses, user axioms, sorry and unchecked evaluation rejected.');

  await checkCalculations({ compile, report, page, reject });
  await checkDerivations({ compile, report, page, reject });

  await compile(arbitrary);
  await context.setOffline(true);
  await page.reload(); await ready(); await compile(arbitrary);
  assert.match(await page.locator('#lean-evaluation').textContent(),/some 7/);
  console.log('Offline reload creates a fresh compiler and checks a new proof.');
  await context.setOffline(false);
  await page.setViewportSize({width:390,height:844});
  assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth > innerWidth),false);
  await page.setViewportSize({width:1400,height:1000});
  console.log('Checking compiler cancellation.');
  await page.locator('#lean-source').fill('import SSC\n#eval IO.sleep 60000');
  await page.locator('#lean-compile').click(); await page.waitForTimeout(500);
  await page.locator('#lean-stop').click();
  await page.waitForFunction(()=>document.getElementById('lean-status').textContent==='Stopped.');
  await ready(); await compile(examples['certified-bounds']);
  assert.deepEqual(requests.filter(r=>r.path.includes('/api/') || /\/study\.(mjs|wasm)/.test(r.path)),[]);
  assert.deepEqual(errors,[]);
  if(process.env.SCREENSHOTS) await page.screenshot({path:resolve(process.env.SCREENSHOTS,'lean-proofs.png'),fullPage:true});
  console.log('Lean proof compiler: generic declarations, statistical proofs, failure modes, offline execution and cancellation passed. No compiler API or numerical backend.');
} catch(error) {
  console.error('Browser status:',await page.locator('#lean-connection-status').textContent());
  console.error('Browser errors:',errors,warnings.slice(-8));
  throw error;
} finally {clearInterval(progressTimer);await context.close();await browser.close();await new Promise(r=>server.close(r));await rm(temporary,{recursive:true,force:true});}
