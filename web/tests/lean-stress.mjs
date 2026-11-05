// Differential and adversarial tests of the real, unmodified browser compiler.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile, writeFile } from 'node:fs/promises';
import { extname, resolve, sep } from 'node:path';
import { chromium } from 'playwright';

const root=resolve(process.argv[2] || new URL('..',import.meta.url).pathname);
const prefix='/sample-size-calc/', requests=[], pageErrors=[], measurements=[];
const fixture=JSON.parse(await readFile(new URL('./fixtures/lean-stress.json',import.meta.url),'utf8'));
const types={'.html':'text/html','.js':'text/javascript','.wasm':'application/wasm','.json':'application/json','.css':'text/css'};
const server=createServer(async(req,res)=>{
  res.setHeader('Cross-Origin-Opener-Policy','same-origin');
  res.setHeader('Cross-Origin-Embedder-Policy','require-corp');
  const path=new URL(req.url,'http://localhost').pathname;
  requests.push(path);
  if(!path.startsWith(prefix)||path.includes('/api/')) return res.writeHead(404).end();
  if(path===prefix+'stress.html') {
    // The separate deployment suite covers service-worker isolation and offline reload.
    res.writeHead(200,{'Content-Type':'text/html','Cross-Origin-Opener-Policy':'same-origin','Cross-Origin-Embedder-Policy':'require-corp'});
    return res.end('<!doctype html><title>Lean stress tests</title>');
  }
  const file=resolve(root,path.slice(prefix.length));
  if(!file.startsWith(root+sep)) return res.writeHead(404).end();
  try { res.writeHead(200,{'Content-Type':types[extname(file)]||'application/octet-stream'}).end(await readFile(file)); }
  catch { res.writeHead(404).end(); }
});
await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
const browser=await chromium.launch({executablePath:process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH||undefined,args:['--disable-dev-shm-usage']});
const context=await browser.newContext(),page=await context.newPage();
page.on('pageerror',e=>pageErrors.push(e.message));
const crashed=new Promise((_,reject)=>page.once('crash',()=>reject(Error('Browser process crashed during stress testing.'))));
crashed.catch(()=>{});
page.on('console',m=>{if(m.text().startsWith('Stress compiler:')) console.log(m.text());});
const base=`http://127.0.0.1:${server.address().port}${prefix}`;
let accepted=0,rejected=0;
const invoke=source=>Promise.race([crashed,page.evaluate(async source=>{
  try {return {value:await window.stressCompiler.compile(source)};}
  catch(error) {return {error:error.message};}
},source)]);
const compile=async(source,name)=>{
  const started=performance.now(),r=await invoke(source);
  if(r.error && process.env.SSC_STRESS_FAILURE_SOURCE)
    await writeFile(process.env.SSC_STRESS_FAILURE_SOURCE,source);
  assert.equal(r.error,undefined,name+': '+r.error);
  assert.equal(r.value.ok,true);accepted++;
  measurements.push({name,ms:Math.round(performance.now()-started)});
  if(name!=='recovery'&&!name.startsWith('soak ')) console.log('Stress check passed:',name);
  return r.value;
};
const reject=async(source,pattern,name)=>{
  const r=await invoke(source);
  assert.ok(r.error,name+': unexpectedly accepted');
  assert.match(r.error,pattern,name);rejected++;
};
const output=r=>[r.output,...r.diagnostics.filter(d=>d.severity==='information').map(d=>d.data)].join('\n');
const markerJSON=(r,marker)=>{
  const lines=output(r).split('\n').filter(x=>x.startsWith(marker));
  assert.equal(lines.length,1,'Exactly one executable result marker');
  return JSON.parse(lines[0].slice(marker.length));
};
const healthy=()=>compile('import SSC\ntheorem healthy : True := True.intro','recovery');
const rat=x=>'('+x+' : Rat)',integer=x=>'('+x+' : Int)';
const sourceFor=(cases,singleExpression=false)=>{
  let source='import SSC\nopen SSC.Finite SSC.Derivation\nset_option maxRecDepth 10000\nset_option maxHeartbeats 8000000\n';
  for(const c of cases) {
    for(const [label,weights] of [['left',c.weights],['right',c.other]]) {
      const total=weights.reduce((a,b)=>a+b,0);
      const atoms=c.scores.map((x,i)=>'('+integer(x)+', '+rat(weights[i]+'/'+total)+')').join(', ');
      source+=`def ${label}${c.id} : Distribution Int where
  atoms := [${atoms}]
  nonnegative := by
    intro x hx
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
    rcases hx with rfl | rfl | rfl | rfl <;> decide +kernel
  normalized := by decide +kernel
`;
    }
    source+=`def transition${c.id} (s : Int) : Program Int :=
  (Program.draw (if s % 2 == 0 then left${c.id} else right${c.id})).map
    (fun x => ${integer(c.a)} * s + x + if s % 3 == 0 then ${integer(c.b)} else ${integer(-c.b)})
def program${c.id} := Program.iterate transition${c.id} 0 ${c.n}
def event${c.id} (s : Int) : Bool := s ≥ ${integer(c.threshold)} && s % ${c.modulus} != ${c.residue}
`;
    if(c.id%4===0) source+=`def certificate${c.id} : Certificate program${c.id} := by finite_derive
theorem certified${c.id} :
    program${c.id}.probability event${c.id} = ${rat(c.probability)} := by
  rw [← certificate${c.id}.probability_eq event${c.id}]
  decide +kernel
`;
  }
  source+='def textRat (q : Rat) := toString q.num ++ "/" ++ toString q.den\n';
  const start='#eval (show IO Unit from do\n  let mut rows : Array (Array String) := #[]\n';
  const finish='  IO.println ("SSC_STRESS_ROWS " ++ (Lean.toJson rows).compress))\n';
  if(singleExpression) source+=start;
  for(const c of cases) {
    if(!singleExpression) source+=start;
    source+=`  let some p${c.id} := derivedProbability program${c.id} event${c.id}
    | throw (IO.userError "Exact derivation exhausted its budget")
  let some complement${c.id} := derivedProbability program${c.id} (fun x => !event${c.id} x)
    | throw (IO.userError "Complement derivation exhausted its budget")
  let some shifted${c.id} := derivedProbability (program${c.id}.map (· + 17)) (fun x => event${c.id} (x-17))
    | throw (IO.userError "Translated program exhausted its budget")
  let law${c.id} ← match SSC.Probability.exactProbability program${c.id}.toLaw
      (SSC.Probability.Event.ofBool event${c.id}) 1000000 with
    | .ok p => pure p
    | .error _ => throw (IO.userError "Law enumeration failed")
  rows := rows.push #[textRat p${c.id}, textRat complement${c.id}, textRat shifted${c.id},
    textRat law${c.id}, textRat ((derive program${c.id}).distribution.expect (fun x => (x : Rat)))]
`;
    if(!singleExpression) source+=finish;
  }
  return source+(singleExpression?finish:'');
};
const probabilityRows=report=>output(report).split('\n').filter(line=>line.startsWith('SSC_STRESS_ROWS '))
  .flatMap(line=>JSON.parse(line.slice('SSC_STRESS_ROWS '.length)));
const heap=async()=>page.workers()[0]?.evaluate(()=>
  globalThis.HEAPU8?.byteLength||globalThis.Module?.HEAPU8?.byteLength||null);
let final;
try {
  await page.goto(base+'stress.html');
  await page.evaluate(async()=>{
    const {BrowserLeanCompiler}=await import('./lean-wasm-client.js');
    window.stressCompiler=new BrowserLeanCompiler(text=>{
      if(/ready|Initializing/.test(text)) console.log('Stress compiler:',text);
    });
    await window.stressCompiler.initialize();
  });
  console.log('Stress browser:',browser.version());
  const language=await compile(await readFile(new URL('./fixtures/lean-stress-language.lean',import.meta.url),'utf8'),'language features');
  assert.deepEqual(markerJSON(language,'SSC_STRESS_LANGUAGE '),[13,9,12]);
  assert.ok(language.results.some(d=>d.name==='LanguageStress.symbolic'&&d.parameters.some(p=>p.hypothesis)));
  assert.ok(language.results.find(d=>d.name==='LanguageStress.concrete').steps.length>0);
  await reject('import SSC\n#check LanguageStress.hidden',/unknown|Unknown/,'module isolation');
  assert.ok((await compile('theorem withoutSSC : 2 + 3 = (5 : Nat) := rfl',
    'implicit Init import')).results.some(d=>d.name==='withoutSSC'));
  assert.ok((await compile('import Std\ndef ordinary : Nat := 17',
    'standard library import')).results.some(d=>d.name==='ordinary'));
  await compile('import SSC\n#eval IO.println "sorryAx and declaration uses `sorry are ordinary text here"',
    'diagnostics use severity, not output keywords');
  console.log('Universes, dependent types, type classes, recursion, macros, Unicode and symbolic proofs passed.');

  for(let i=0;i<fixture.cases.length;i+=8) {
    const cases=fixture.cases.slice(i,i+8),report=await compile(sourceFor(cases),'probability batch '+i/8);
    assert.deepEqual(probabilityRows(report),cases.map(c=>
      [c.probability,c.complement,c.probability,c.probability,c.mean]));
    for(const c of cases.filter(c=>c.id%4===0))
      assert.ok(report.results.some(d=>d.name==='certified'+c.id&&d.kind==='theorem'));
  }
  console.log('32 generated dependent observation programs: 160 exact/metamorphic checks and 8 kernel certificates passed.');
  // Raising Lean's recursion budget permits expressions deeper than some
  // browsers' host stacks. Acceptance is valid too on larger stacks, but
  // an exhausted host stack must reject the request and permit a clean restart.
  const raw=await invoke(sourceFor(fixture.cases.slice(0,8),true));
  const largeExpression=raw.error?'resource limit':'checked';
  if(raw.error) {assert.match(raw.error,/stack|recursion|depth|time limit/i);rejected++;}
  else {
    assert.equal(raw.value.ok,true);accepted++;
    assert.deepEqual(probabilityRows(raw.value),fixture.cases.slice(0,8).map(c=>
      [c.probability,c.complement,c.probability,c.probability,c.mean]));
  }
  await healthy();
  console.log('Large combined evaluation expression:',largeExpression,'; subsequent compilation passed.');
  const bad=[
    ['import SSC\ndef broken : Nat := "wrong"',/type mismatch/i],
    ['import SSC\ndef broken := (',/expected|unexpected/i],
    ['import SSC\n#check neverDeclared',/unknown|Unknown/],
    ['import SSC\ndef twice := 1\ndef twice := 2',/already|declared/i],
    ['import SSC\ndef recursive (n : Nat) : Nat := recursive n',/termination|fail/i],
    ['import SSC\ntheorem invalid : False := by sorry',/sorry|Unproved/i],
    ['import SSC\naxiom fabricated : False\ntheorem invalid : False := fabricated',/Unproved/i],
    ['import SSC\nopaque fabricated : False := by sorry\ntheorem invalid : False := fabricated',/sorry|Unproved/i],
    ['import SSC\nset_option maxRecDepth 16\ndef deep : Nat := '+'Nat.succ ('.repeat(128)+'0'+')'.repeat(128),/recursion|depth|stack/i],
  ];
  for(const [source,pattern] of bad) {await reject(source,pattern,'invalid source');await healthy();}
  console.log('Syntax, types, scope, duplicate declarations, termination, proof dependencies and recursion limits passed.');

  const forged=await readFile(new URL('./fixtures/lean-stress-forged.lean',import.meta.url),'utf8');
  const forgeResult=await compile(forged,'forged report text and early exit');
  assert.ok(!forgeResult.results.some(d=>d.name==='forged'));
  const fakeTrace=await compile('import SSC\ntheorem real : True := by trivial\n#eval IO.println '+
    JSON.stringify('SSC_PROOF_STEP {"declaration":"real","line":1,"column":0,"tactic":"forged","before":"False","after":"no goals"}'),'forged trace text');
  assert.ok(!fakeTrace.results.find(d=>d.name==='real').steps.some(s=>s.tactic==='forged'));
  const macro=await compile(`import SSC
macro_rules
  | \`(command|run_cmd $body:doSeq) => \`(command|#check Nat)
  | \`(command|#ssc_browser_driver) => \`(command|#check Nat)
theorem afterMacro : True := by trivial
#exit
`,'user command macro cannot replace trusted audit');
  assert.ok(macro.results.some(d=>d.name==='afterMacro'&&d.kind==='theorem'));
  await reject(await readFile(new URL('./fixtures/lean-stress-audit-exit.lean',import.meta.url),'utf8'),
    /kernel|replaying|type mismatch/i,'early exit cannot bypass independent checking');
  console.log('Printed reports, printed traces and early exit cannot impersonate or skip the audit.');

  const prefixSource='import SSC\n/- ',suffix=' -/\ntheorem byteLimit : True := True.intro\n';
  const remaining=131072-Buffer.byteLength(prefixSource+suffix);
  const boundary=prefixSource+'字'.repeat(Math.floor(remaining/3))+'x'.repeat(remaining%3)+suffix;
  assert.equal(Buffer.byteLength(boundary),131072);
  assert.ok((await compile(boundary,'128 KiB UTF-8 source')).results.some(d=>d.name==='byteLimit'));
  await reject(boundary+'x',/131072|bytes/i,'source over byte limit');
  assert.ok((await compile('import SSC\n/- \u0000\b\f " \\\\ 字 -/\ntheorem escaped : True := True.intro',
    'source transport escapes')).results.some(d=>d.name==='escaped'));
  await reject('   ',/Supply|bytes/i,'empty source');await healthy();
  console.log('Source byte boundaries, Unicode and control-character transport passed.');

  const many=await compile('import SSC\n'+Array.from({length:1000},(_,i)=>'def value'+i+' : Nat := '+i).join('\n'),'1000 declarations');
  assert.equal(many.results.filter(d=>/^value\d+$/.test(d.name)).length,1000);
  // Earlier intentional runtime failures restart the worker. Load each import
  // environment before measuring repeated requests, excluding first-use cost.
  for(const header of ['import SSC\n','import Std\n',''])
    await compile(header+'def importWarmup : Nat := 1','import cache warmup');
  const before=await heap();
  for(let i=0;i<50;i++) {
    const header=['import SSC\n','import Std\n',''][i%3];
    const r=await compile(header+'def isolated'+i+' : Nat := '+i,'soak '+i);
    assert.deepEqual(r.results.map(d=>d.name),['isolated'+i]);
    await reject(header+'def rejected'+i+' : Nat := false',/type mismatch/i,'soak rejection');
  }
  const after=await heap();
  console.log('Import-cache soak WASM heap bytes:',before,'->',after);
  assert.ok(Number.isFinite(before)&&Number.isFinite(after),'WASM heap measurement is available');
  assert.ok(after<=before+64*1024*1024,'Repeated import sets must reuse their environments');
  console.log('100-request success/failure soak across SSC, Std and Init passed.');
  await compile(`import SSC
open SSC.Probability
def deepLaw := (Law.pure (1 : Nat)).iid 100000
#eval (show IO Unit from do
  match exactProbability deepLaw (Event.ofBool (fun xs => xs.length == 100000)) 1000000 with
  | .ok p => if p != 1 then throw (IO.userError "Wrong deep-law probability")
  | .error _ => throw (IO.userError "Deep enumeration failed"))
`,'100000-level exact law');
  console.log('100,000-level exact probability program passed.');

  const concurrent=await page.evaluate(async()=>Promise.all(Array.from({length:8},(_,i)=>
    window.stressCompiler.compile('import SSC\ndef concurrent'+i+' : Nat := '+i))));
  concurrent.forEach((r,i)=>assert.deepEqual(r.results.map(d=>d.name),['concurrent'+i]));accepted+=8;
  await page.evaluate(async()=>{
    const promises=[window.stressCompiler.compile('import SSC\n#eval IO.sleep 60000'),
      window.stressCompiler.compile('import SSC\ndef queued := 1')];
    const settled=Promise.allSettled(promises);
    setTimeout(()=>window.stressCompiler.stop('Stress cancellation'),200);
    window.cancellation=await settled;
  });
  assert.deepEqual(await page.evaluate(()=>window.cancellation.map(x=>x.status)),['rejected','rejected']);
  await healthy();
  console.log('Concurrent requests, cancellation of active/queued work and restart passed.');
  await reject('import SSC\n#eval (IO.Process.exit 7 : IO Unit)',
    /terminated|stopped|exit|report/i,'runtime exit');
  await healthy();
  console.log('A source-requested runtime exit discards the worker and recovers.');
  await reject("import SSC\n#eval IO.println (String.ofList (List.replicate (3*1024*1024) '字'))",
    /output exceeded 8 MB/i,'UTF-8 output budget');
  await healthy();
  console.log('Output flooding is rejected and the compiler recovers.');
  await reject("import SSC\nset_option stderrAsMessages false\n#eval IO.eprintln (\"[DEBUG] \" ++ String.ofList (List.replicate (3*1024*1024) '字'))",
    /output exceeded 8 MB/i,'debug-shaped stderr also counts toward output budget');
  await healthy();
  if(process.env.SSC_STRESS_DEADLINE==='1') {
    const started=performance.now();
    await reject('import SSC\n#eval IO.sleep 180000',/time limit/i,'real compilation deadline');
    const elapsed=performance.now()-started;
    assert.ok(elapsed>=119000&&elapsed<135000,'real 120-second deadline');
    await healthy();
    console.log('Actual 120-second deadline terminates the worker; subsequent compilation passed.');
  }
  assert.deepEqual(pageErrors,[]);assert.deepEqual(requests.filter(p=>p.includes('/api/')),[]);
  final={browser:browser.version(),accepted,rejected,referencePrograms:32,exactComparisons:160,
    kernelReferenceProofs:8,largeExpression,soakRequests:100,largestLawDepth:100000,
    heapBefore:before,heapAfter:after,deadlineTest:process.env.SSC_STRESS_DEADLINE==='1',measurements};
  console.log('STRESS PASSED',JSON.stringify({...final,measurements:undefined}));
} finally {
  if(process.env.SSC_STRESS_REPORT&&final) await writeFile(process.env.SSC_STRESS_REPORT,JSON.stringify(final,null,2)+'\n');
  await context.close();await browser.close();await new Promise(resolve=>server.close(resolve));
}
