import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { execFileSync } from 'node:child_process';
import { examples } from '../lean-examples.js';

const executable = process.env.SSC_DESIGN_NATIVE;
const M = executable ? null : await (await import('../engine/ssc_web.mjs')).default();
function design(request) {
  if (executable) {
    let output;
    try { output = execFileSync(executable, ['design', '-'], { input: JSON.stringify(request), encoding: 'utf8' }); }
    catch (error) { output = error.stdout; }
    const answer = JSON.parse(output);
    if (!answer.ok) throw new Error(answer.error);
    return answer.result;
  }
  const text = JSON.stringify(request), p = M._malloc(M.lengthBytesUTF8(text) + 1);
  assert.ok(p); let out;
  try {
    M.stringToUTF8(text, p, M.lengthBytesUTF8(text) + 1);
    out = M._ssc_design_json(p, 0); assert.ok(out);
    const answer = JSON.parse(M.UTF8ToString(out));
    if (!answer.ok) throw new Error(answer.error);
    return answer.result;
  } finally { if (out) M._ssc_free(out); M._free(p); }
}

test('editor examples are the actual Lean files', async () => {
  for (const name of ['normal', 'influence']) assert.equal(examples[name], await readFile(new URL(`../../examples/design/${name}.lean`, import.meta.url), 'utf8'));
});
test('Lean declarations compile in C++ and solve all three quantities', () => {
  for (const name of ['normal', 'influence', 'cluster']) {
    const count = name === 'influence' ? 126 : 63;
    for (const mode of ['sample_size', 'power', 'effect']) {
      const request = { op: 'compile', source: examples[name], solveFor: mode };
      if (mode !== 'sample_size') request.n = count;
      const c = design(request), r = design(c.plan);
      assert.equal(c.leanKernelChecked, false);
      assert.match(c.cpp, /int main\(\)/);
      if (mode === 'sample_size') assert.equal(r.selectedN, count);
      if (mode === 'power') assert.ok(Math.abs(r.power - .8013023941055797) < 2e-13);
      if (mode === 'effect') { assert.ok(r.detectableMagnitude < .5); assert.ok(r.powerAtUpper >= .8); }
    }
  }
});
test('compiler rejects unsupported Lean and incompatible study fields', () => {
  for (const source of [examples.normal + '\n#eval IO.println "test"', examples.normal.replace('sd := 1', 'sd := 0'), examples.normal.replace('sd := 1', 'variance := 1')])
    assert.throws(() => design({ op: 'compile', source }));
});
test('third tab has separate navigation and offline files', async () => {
  for (const page of ['index.html', 'dag.html', 'lean.html']) assert.match(await readFile(new URL(`../${page}`, import.meta.url), 'utf8'), /href="lean\.html\?v=dev"/);
  const sw = await readFile(new URL('../sw.js', import.meta.url), 'utf8');
  for (const file of ['lean.html', 'lean.css', 'lean.js', 'lean-examples.js']) assert.ok(sw.includes(`${file}?v=dev`));
});
