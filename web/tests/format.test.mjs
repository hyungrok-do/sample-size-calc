// The page formats numbers as the command line does (python/ssc/formatting.py): the same values go through both
// and must come out the same.  A negative zero is the one difference by design: the page shows it as 0.
import assert from 'node:assert/strict';
import test from 'node:test';
import { spawnSync } from 'node:child_process';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { cell } from '../format.js';

const repo = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const python = process.env.PYTHON ?? (process.platform === 'win32' ? 'python' : 'python3');
const FORMATS = ['num', 'prob', 'int', 'pct', 'alpha', 'ratio'];

function values() {
  const vals = [0, 1, -1, 0.5, 2.5, -2.5, 1.5, 0.03125, -0.03125, 0.00005, 0.0001, 0.00009999, 0.000099996, 0.00012345,
    1e-5, 1.2345e-7, -3.2e-5, 999999.5, 1e6, 1234567.891, -1234567.5, 1e15, 1e16, 123456789012345678, 1e21, 1e22, 12345,
    9999, 10000, -10000, 99999.5, 0.99995, 0.999949999, 0.12345, 0.123449999, 1.00005, 2.675, 1.0000000001, -0.00004,
    -1e-10, 5e-324, 1.7976931348623157e308, 0.125, 0.375, 0.0625, 0.00625, 1 / 3, 2 / 3, 33.3333, 0.84375, 0.15625,
    100.5, 0.4, -0.4, -0.5, 1e-4 - 1e-20, 0.29, 0.2949999, 12.5, 0.000125, 3.0000001, 3663969941044.34375,
    2 ** 52 + 0.5, 2 ** 51 + 0.25, 123456.03125, -123456.03125, null, true, false, 'text', NaN, Infinity, -Infinity];
  let x = 12345;  // xorshift: the same values on every run
  const rnd = () => { x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0; return x / 4294967296; };
  for (let i = 0; i < 4000; i++) {
    const e = Math.floor(rnd() * 24) - 12;
    vals.push((rnd() < 0.2 ? -1 : 1) * rnd() * 10 ** e);
    vals.push(Math.round(rnd() * 10 ** (e + 5)) / 2 ** Math.floor(rnd() * 8));
  }
  // exact ties at every number of places, above 2^53 too
  for (let i = 0; i < 2000; i++) vals.push((Math.floor(rnd() * 2 ** 40) + 0.5) / 2 ** Math.floor(rnd() * 12));
  return vals;
}

test('numbers look the same on the page and on the command line', () => {
  const vals = values();
  // JSON has no NaN or infinities: those travel as {"float": "..."}
  const sent = vals.map((v) => (typeof v === 'number' && !Number.isFinite(v) ? { float: String(v) } : v));
  const code = [
    'import json, sys',
    'from ssc.formatting import fmt_cell',
    'vals, fmts = json.load(sys.stdin)',
    'vals = [float(v["float"]) if isinstance(v, dict) else v for v in vals]',
    'print(json.dumps([[fmt_cell(v, f) for f in fmts] for v in vals]))',
  ].join('\n');
  const run = spawnSync(python, ['-c', code], {
    input: JSON.stringify([sent, FORMATS]),
    env: { ...process.env, PYTHONPATH: resolve(repo, 'python') },
    encoding: 'utf8',
    maxBuffer: 64 * 1024 * 1024,
  });
  assert.equal(run.status, 0, run.stderr);
  const expected = JSON.parse(run.stdout);
  const differ = [];
  vals.forEach((v, i) => FORMATS.forEach((f, j) => {
    const page = cell(v, f);
    if (page !== expected[i][j]) differ.push(`${String(v)} as ${f}: page ${page}, command line ${expected[i][j]}`);
  }));
  assert.deepEqual(differ.slice(0, 10), [], `${differ.length} of ${vals.length * FORMATS.length} differ`);
  assert.equal(cell(-0, 'num'), '0');
});
