// Numbers, as the command line shows them (python/ssc/formatting.py): the same rounding and exponents.
// web/tests/format.test.mjs compares them with the command line's over many values.
export const trimZeros = (s) => (s.includes('.') ? s.replace(/0+$/, '').replace(/\.$/, '') : s);
export const grouped = (s) => s.replace(/\B(?=(\d{3})+(?!\d))/g, ',');

// toFixed as Python formats: an exact tie goes to the even digit (0.03125 is exactly 1/32, and to 4 places is
// 0.0312), while 0.00625, a little above a tie as a double, rounds up as toFixed does; and no exponent above 1e21
export function fixed(v, d) {
  const a = Math.abs(v);
  if (a >= 1e21) return (v < 0 ? '-' : '') + BigInt(a).toString() + (d ? `.${'0'.repeat(d)}` : '');
  const places = Math.min(100, d + 30);  // the double's exact decimal digits (30 more settle any tie)
  const [whole, exact = ''] = a.toFixed(places).split('.');
  if (exact[d] === '5' && /^0*$/.test(exact.slice(d + 1))) {
    // the digits kept, plus one in the last place when it is odd (in whole numbers: above 2^53 too)
    const kept = BigInt(whole + exact.slice(0, d));
    const digits = (kept % 2n ? kept + 1n : kept).toString().padStart(d + 1, '0');
    return (v < 0 ? '-' : '') + (d ? `${digits.slice(0, -d)}.${digits.slice(-d)}` : digits);
  }
  return v.toFixed(d);
}

// Python's '{:.4g}': fixed notation for exponents from -4 to 3, otherwise scientific with a two-digit exponent
export function general(v, p = 4) {
  if (v === 0) return '0';
  const [mantissa, e] = v.toExponential(p - 1).split('e');
  const exp = Number(e);
  if (exp >= -4 && exp < p) return trimZeros(fixed(v, p - 1 - exp));
  return `${trimZeros(mantissa)}e${exp < 0 ? '-' : '+'}${String(Math.abs(exp)).padStart(2, '0')}`;
}

export function smart(v, maxDec = 4) {
  if (!Number.isFinite(v)) return '—';
  if (Number.isInteger(v) && Math.abs(v) < 1e15) return Math.abs(v) >= 10000 ? grouped(fixed(v, 0)) : fixed(v, 0);
  const a = Math.abs(v);
  if (a >= 1e6 || a < 1e-4) return general(v, 4);
  const s = trimZeros(fixed(v, maxDec));
  return s === '-0' || s === '' ? '0' : s;
}

export function cell(v, fmt) {
  if (v == null) return '—';
  if (typeof v === 'string') return v;
  if (typeof v === 'boolean') return v ? 'Yes' : 'No';
  if (!Number.isFinite(v)) return '—';
  if (fmt === 'prob') return fixed(v, 4);
  if (fmt === 'int') return grouped(fixed(v, 0));
  if (fmt === 'pct') return `${smart(100 * v, 2)}%`;
  return smart(v);
}
