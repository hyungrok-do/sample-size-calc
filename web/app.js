// The calculator in a web page, on its own engine: pick a method, fill its inputs (a form generated from the
// method's schema) and calculate; and a speed test of the engine in this browser.
import { Engine } from './engine.js?v=dev';
import { mountGuide } from './guide.js?v=dev';
import { visible, prepareSelection, selectionHash, parseSelectionHash } from './selection.js?v=dev';
import { cell } from './format.js?v=dev';

const $ = (id) => document.getElementById(id);
const engine = new Engine();
const state = { catalog: null, proc: null, solve: '', fields: new Map(), busy: false, form: 0 };
const schemas = new Map();  // procedure id → its description, so opening a method never waits on a calculation
const MAX_ROWS = 1000;      // rows of an answer drawn in the table
let guide;
let mode = 'browse';
let guidedForm = false;

function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, v] of Object.entries(attrs)) {
    if (v === false || v == null) continue;
    if (key === 'class') node.className = v;
    else if (key.startsWith('on')) node.addEventListener(key.slice(2), v);
    else node.setAttribute(key, v === true ? '' : v);
  }
  node.append(...children.flat().filter((c) => c != null && c !== false));
  return node;
}

// Times and sizes on the page (numbers in results are formatted by format.js, as the command line does)
const ms =(t) => (t < 10 ? `${t.toFixed(1)} ms` : t < 1000 ? `${Math.round(t)} ms` : `${(t / 1000).toFixed(1)} s`);
const mb = (bytes) => `${(bytes / 1e6).toFixed(bytes < 1e6 ? 2 : 1)} MB`;

// ---------------------------------------------------------------------------------------------------------
// The engine's status line
// ---------------------------------------------------------------------------------------------------------
function describe(info) {
  const parts = [`Engine ${info.version}`, info.threads > 1 ? `${info.threads} threads` : '1 thread'];
  const d = info.download;
  if (d && d.encoded) {
    // a revalidated or cached file arrives in a few hundred bytes of headers, or none
    if (d.transfer < d.encoded / 2) parts.push('from the cache');
    else parts.push(`${mb(d.encoded)} download${d.decoded > d.encoded * 1.05 ? ` (${mb(d.decoded)} unpacked)` : ''}`);
  }
  parts.push(`ready in ${ms(info.loadMs)}`);
  return parts.join(' · ');
}

async function showStatus() {
  try {
    $('engine-status').textContent = describe(await engine.ready);
    return true;
  } catch (e) {
    $('engine-status').textContent = state.catalog
      ? `The engine could not restart: ${e.message}. It tries again when the connection returns; reloading the page also works.`
      : `The engine did not load: ${e.message}`;
    return false;
  }
}

// ---------------------------------------------------------------------------------------------------------
// Methods: the catalog by category, or what a search finds
// ---------------------------------------------------------------------------------------------------------
function methodButton(p) {
  const current = state.proc?.id === p.id;
  return el('li', {}, el('button', {
    type: 'button', class: 'method', 'data-id': p.id, 'aria-current': current ? 'true' : null, onclick: () => choose(p.id),
  }, p.name));
}

function renderList(hits, note) {
  const nav = $('method-list');
  nav.replaceChildren();
  if (hits) {
    if (note) nav.append(el('p', { class: 'muted small' }, note));
    if (!hits.length) nav.append(el('p', { class: 'muted small' }, 'No method matches these words.'));
    nav.append(el('ul', {}, hits.map(methodButton)));
    return;
  }
  for (const category of state.catalog.categories) {
    const procs = state.catalog.procedures.filter((p) => p.category === category);
    if (!procs.length) continue;
    nav.append(el('h3', {}, category), el('ul', {}, procs.map(methodButton)));
  }
}

let searchTimer = 0;
let searchSeq = 0;
let searchedWhileBusy = false;
function onSearch() {
  clearTimeout(searchTimer);
  searchTimer = setTimeout(search, 120);
}

async function search() {
  if (!state.catalog) return;  // still loading: the list is drawn with this query once the catalog is in
  const query = $('search').value.trim();
  const seq = ++searchSeq;
  if (!query) {
    renderList(null);
    return;
  }
  if (state.busy) {  // the engine is calculating: match names and keywords here, and search again when it is done
    searchedWhileBusy = true;
    renderList(quickSearch(query), 'Matching names and keywords while the calculation runs.');
    return;
  }
  try {
    const { value } = await engine.call('search', { query, limit: 20 });
    if (seq === searchSeq) renderList(value.hits, value.corrected ? `Showing results for “${value.corrected}”` : '');
  } catch {
    // stopped: the next keystroke searches again
  }
}

const plainWords = (text) => text.toLowerCase().replace(/[^\p{L}\p{N}\s]/gu, '');
function quickSearch(query) {
  const words = plainWords(query).split(/\s+/).filter(Boolean);
  return state.catalog.procedures.filter((p) => {
    const text = plainWords([p.name, p.id, p.category, p.subcategory, ...(p.keywords ?? [])].join(' '));
    return words.every((w) => text.includes(w));
  }).slice(0, 20);
}

// ---------------------------------------------------------------------------------------------------------
// One method: solve-for choices and the generated form
// ---------------------------------------------------------------------------------------------------------
function value(id) {
  const f = state.fields.get(id);
  if (!f) return '';
  if (f.p.type === 'bool') return f.input.checked ? 'true' : 'false';
  if (f.p.type === 'select' && f.options === null) return f.initial;
  return f.input.value;
}

const lookup = (key) => (key === '$solve' ? state.solve : value(key));

// Shows the inputs that apply to the current choices and offers the options that do.
function refresh() {
  for (let pass = 0; pass < 3; pass++) {
    let changed = false;
    for (const f of state.fields.values()) {
      if (f.p.type !== 'select') continue;
      const wanted = f.p.options.filter((o) => visible(o.visible, lookup));
      const keys = wanted.map((o) => o.value).join('\n');
      if (keys === f.options) continue;
      const was = value(f.p.id);
      f.input.replaceChildren(...wanted.map((o) => el('option', { value: o.value }, o.label)));
      const keep = wanted.find((o) => o.value === was) ?? wanted.find((o) => o.value === f.p.default) ?? wanted[0];
      f.input.value = keep ? keep.value : '';
      f.options = keys;
      changed = true;
    }
    if (!changed) break;
  }
  for (const f of state.fields.values()) f.wrap.hidden = !visible(f.p.visible, lookup);
  for (const fs of $('fields').querySelectorAll('fieldset')) {
    fs.hidden = ![...fs.querySelectorAll('.field')].some((w) => !w.hidden);
  }
}

// The phone keyboard for a field: digits only for a single whole number that cannot be negative; the full
// keyboard otherwise (lists and ranges need spaces and words, and iOS's decimal pad has no minus sign)
function inputmode(p) {
  if (p.type !== 'number' || p.multi) return null;
  return p.integer && p.min != null && p.min >= 0 ? 'numeric' : null;
}

function field(p, initial = p.default ?? '') {
  const id = `in-${p.id}`;
  let input;
  if (p.type === 'select') input = el('select', { id, onchange: refresh });
  else if (p.type === 'bool') input = el('input', { id, type: 'checkbox', checked: initial === 'true', onchange: refresh });
  else {
    input = el('input', {
      id, type: 'text', value: initial, autocomplete: 'off', spellcheck: 'false', inputmode: inputmode(p),
      oninput: () => input.removeAttribute('aria-invalid'),
    });
  }
  // two lines of help, all of it on a click (and on hover, and for screen readers): it never grows by itself
  // under the pointer
  const help = p.help
    ? el('p', { class: 'help', id: `${id}-help`, title: p.help, onclick: (e) => e.currentTarget.classList.toggle('open') }, p.help)
    : null;
  if (help) input.setAttribute('aria-describedby', help.id);
  const label = el('label', { for: id }, p.label);
  const wrap = el('div', { class: 'field' }, p.type === 'bool' ? el('div', { class: 'check' }, input, label) : [label, input], help);
  return { p, wrap, input, options: null, initial };
}

function buildForm(values = {}) {
  const box = $('fields');
  // Chromium's style engine can keep a reference to an option of a removed <select>, and through it the
  // whole old form: emptying the selects first leaves it at most an option
  for (const f of state.fields.values()) if (f.p.type === 'select') f.input.replaceChildren();
  box.replaceChildren();
  state.fields.clear();
  state.form++;
  for (const group of state.proc.groups) {
    const params = state.proc.params.filter((q) => q.group === group);
    if (!params.length) continue;
    const grid = el('div', { class: 'grid' });
    for (const p of params) {
      const f = field(p, values[p.id]);
      state.fields.set(p.id, f);
      grid.append(f.wrap);
    }
    box.append(el('fieldset', {}, el('legend', {}, group), grid));
  }
  refresh();
}

function renderSolve() {
  const box = $('solve');
  box.replaceChildren(el('span', { class: 'solve-label', id: 'solve-label' }, 'Solve for'));
  for (const s of state.proc.solve_for) {
    box.append(el('label', { class: 'seg' },
      el('input', {
        type: 'radio', name: 'solve', value: s.value, checked: s.value === state.solve,
        onchange: () => {
          state.solve = s.value;
          state.form++;
          refresh();
          clearResult();
        },
      }),
      el('span', {}, s.label)));
  }
}

async function schema(id) {
  if (!schemas.has(id)) schemas.set(id, (await engine.call('procedure', id)).value);
  return schemas.get(id);
}

// The address names the method on screen: a method chosen from the list is a step in the browser's history.  It
// changes only when it must, and a refusal (Safari allows 100 changes in 10 seconds) leaves the address as it was
// rather than failing the action.
function address(url, push = false) {
  if (url === location.pathname + location.search + location.hash) return;
  try { history[push ? 'pushState' : 'replaceState'](null, '', url); } catch { /* the address stays as it was */ }
}
const known = (id) => Boolean(state.catalog?.procedures.some((p) => p.id === id));
function choose(id) {
  // the method already open keeps what was typed into it
  if (id === state.proc?.id && !$('method').hidden) { setMode('browse'); return; }
  navigate({ procedure: id, inputs: {} });
}

function setMode(next) {
  mode = next;
  $('workspace').dataset.mode = mode;
  $('browse-methods').hidden = mode !== 'browse';
  $('guide').hidden = mode !== 'guide';
  $('mode-guide').setAttribute('aria-pressed', String(mode === 'guide'));
  $('mode-browse').setAttribute('aria-pressed', String(mode === 'browse'));
  $('method').hidden = !state.proc || (mode === 'guide' && !guidedForm);
  if (mode === 'guide') $('navigation-error').hidden = true;  // it was about a method, not the guide
}

async function navigate(selection, fromGuide = false, push = true) {
  $('navigation-error').hidden = true;
  try {
    if (!(await openMethod(selection, fromGuide))) return false;
    if (push) address(location.pathname + location.search + selectionHash(selection), true);
    setMode(fromGuide ? 'guide' : 'browse');
    $('m-name').focus();
    return true;
  } catch (error) {
    if (fromGuide) throw error;
    state.proc = null;
    guidedForm = false;
    guide?.clearApplied();
    setMode('browse');
    $('method').hidden = true;
    $('navigation-error').textContent = error.message;
    $('navigation-error').hidden = false;
    return false;
  }
}

let openSeq = 0;
async function openMethod(selection, fromGuide = false) {
  const seq = ++openSeq;
  state.form++;  // invalidate an answer before waiting for a new method
  clearResult();
  let proc;
  try {
    proc = await schema(selection.procedure);
  } catch (e) {
    if (!e.stopped && seq === openSeq) throw e;
    return false;
  }
  if (seq !== openSeq) return false;  // a newer choice is already on its way
  const prepared = prepareSelection(proc, selection);
  // Preserve numerical inputs only when their context (design, effect scale, distribution, etc.) is the
  // same. Test/power-calculation choices can change without reinterpreting these inputs. In particular,
  // do not carry a risk-difference margin into a risk-ratio design, or a within-person CV into parallel BE.
  const sameContext = state.proc?.id === proc.id && proc.params
    .filter((p) => p.type === 'select' && !['test', 'method'].includes(p.id))
    .every((p) => value(p.id) === prepared.values[p.id]);
  if (fromGuide && sameContext) {
    for (const f of state.fields.values()) {
      if (f.p.type !== 'select' && f.p.type !== 'bool') prepared.values[f.p.id] = value(f.p.id);
    }
  }
  state.proc = proc;
  state.solve = prepared.solve;
  $('m-category').textContent = [proc.category, proc.subcategory].filter(Boolean).join(' · ');
  $('m-name').textContent = proc.name;
  $('m-summary').textContent = proc.summary;
  $('m-sim').hidden = !proc.simulation;
  renderSolve();
  buildForm(prepared.values);
  for (const [key, wanted] of Object.entries(selection.inputs ?? {})) {
    if (value(key) !== wanted || state.fields.get(key).wrap.hidden) {
      $('method').hidden = true;
      throw new Error('The selected analysis settings could not be applied. Please review the study choices.');
    }
  }
  clearResult();
  guidedForm = fromGuide;
  if (!fromGuide) guide?.clearApplied();
  $('m-source').hidden = !fromGuide;
  $('m-source').textContent = 'From your study choices. Review the example values below and replace them with your study inputs.';
  $('method').hidden = false;
  for (const b of document.querySelectorAll('#method-list .method')) {
    if (b.dataset.id === proc.id) b.setAttribute('aria-current', 'true');
    else b.removeAttribute('aria-current');
  }
  return true;
}

// ---------------------------------------------------------------------------------------------------------
// Calculating, stopping and the answer
// ---------------------------------------------------------------------------------------------------------
// `what` replaces the line next to Calculate; undefined leaves it (the speed test keeps the answer's timing).
function setBusy(busy, what) {
  state.busy = busy;
  // the buttons stay enabled (a disabled button drops the keyboard focus) and say that they wait
  for (const id of ['calc', 'bench']) $(id).setAttribute('aria-disabled', String(busy));
  $('stop').hidden = !busy;
  $('activity-stop').hidden = !busy;
  if (what !== undefined) $('timing').textContent = what;
  if (!busy && searchedWhileBusy) {
    searchedWhileBusy = false;
    search();
  }
}

let drawing = 0;  // the answer whose rows are being drawn: a newer one stops it

function clearResult() {
  drawing++;
  $('result').replaceChildren();
  $('announce').textContent = '';
  for (const f of state.fields.values()) f.input.removeAttribute('aria-invalid');
  if (!state.busy) $('timing').textContent = '';
}

function showError(message) {
  drawing++;
  $('result').replaceChildren(el('div', { class: 'error', role: 'alert' }, el('p', {}, message)));
}

function labelOf(id) {
  return state.proc.params.find((p) => p.id === id)?.label ?? id;
}

// A table row of an answer; its cells are an object (key → value) or an array in column order.
function rowOf(r, cols) {
  if (r._error) return el('tr', {}, el('td', { colspan: cols.length, class: 'row-error' }, `Not computed: ${r._error}`));
  return el('tr', {}, cols.map((c, i) => el('td', {}, cell(Array.isArray(r) ? r[i] : r[c.key], c.fmt))));
}

// Rows a hundred per frame, so a long table never holds the page still.
function drawRows(tbody, rows, cols) {
  const job = drawing;
  let i = 0;
  const step = () => {
    if (job !== drawing) return;
    const batch = document.createDocumentFragment();
    for (const end = Math.min(rows.length, i + 100); i < end; i++) batch.append(rowOf(rows[i], cols));
    tbody.append(batch);
    if (i < rows.length) requestAnimationFrame(step);
  };
  step();
}

function table(cols, rows, label) {
  const head = el('tr', {}, cols.map((c) => el('th', { scope: 'col', title: c.tip ?? null }, c.label)));
  const body = el('tbody');
  const wrap = el('div', { class: 'table-wrap', tabindex: '0', role: 'region', 'aria-label': label },
    el('table', {}, el('thead', {}, head), body));
  drawRows(body, rows, cols);
  return wrap;
}

// The extra tables of one scenario: group-sequential boundaries, Simon's designs, the power of each arm ...
function details(rows) {
  const which = rows.map((r, i) => i).filter((i) => rows[i]._tables?.length);
  if (!which.length) return [];
  const box = el('div', { class: 'details-tables' });
  const show = (i) => box.replaceChildren(...rows[i]._tables.flatMap((t) => [el('h3', {}, t.title), table(t.columns, t.rows, t.title)]));
  const parts = [];
  if (rows.length > 1) {
    const pick = el('select', { id: 'details-row', onchange: (e) => show(Number(e.target.value)) },
      which.map((i) => el('option', { value: String(i) }, `Scenario ${i + 1}`)));
    parts.push(el('div', { class: 'details' }, el('label', { for: 'details-row' }, 'Details for'), pick));
  }
  show(which[0]);
  return [...parts, box];
}

// A message as shown: the engine quotes typed text briefly, and this keeps any other long word from costing the
// page seconds of layout in an error box (which breaks words anywhere).
const brief = (text) => (text.length > 500 ? `${text.slice(0, 500)}…` : text);

function renderResult(res, engineMs) {
  const out = $('result');
  drawing++;
  out.replaceChildren();
  $('timing').textContent = `${ms(engineMs)} in the engine${res.threads > 1 ? ` on ${res.threads} threads` : ''}`;
  if (!res.ok) {
    const errors = Object.entries(res.input_errors ?? {}).map(([id, m]) => {
      state.fields.get(id)?.input.setAttribute('aria-invalid', 'true');
      return el('li', {}, el('strong', {}, `${labelOf(id)}: `), brief(String(m)));
    });
    out.append(el('div', { class: 'error', role: 'alert' },
      el('p', {}, brief(res.error || 'The inputs need attention.')), errors.length ? el('ul', {}, errors) : null));
    return;
  }
  const rows = res.rows ?? [];
  const total = res.rows_total ?? rows.length;
  const statement = rows.map((r) => r._statement).find(Boolean);
  if (statement) out.append(el('p', { class: 'statement' }, statement));
  $('announce').textContent = statement ?? (rows[0]?._error ? `Not computed: ${rows[0]._error}` : 'The answer is ready.');
  out.append(table(res.columns ?? [], rows, 'Results'));
  if (total > rows.length) {
    out.append(el('p', { class: 'muted small' }, `Showing the first ${rows.length.toLocaleString('en-US')} of ${total.toLocaleString('en-US')} scenarios, one per row.`));
  } else if (rows.length > 1) {
    out.append(el('p', { class: 'muted small' }, `${rows.length} scenarios, one per row.`));
  }
  const dropout = (res.dropout?.rows ?? []).filter((r) => r.dropout > 0);
  if (dropout.length) {
    out.append(el('h3', {}, 'Sample size allowing for dropout'), table(res.dropout.columns, dropout, 'Sample size allowing for dropout'));
  }
  out.append(...details(rows));
  if (res.notes?.length) out.append(el('ul', { class: 'notes' }, res.notes.map((n) => el('li', {}, n))));
}

async function calculate(event) {
  event?.preventDefault();
  if (state.busy || !state.proc) return;
  const inputs = {};
  for (const f of state.fields.values()) if (!f.wrap.hidden) inputs[f.p.id] = value(f.p.id);
  const request = { procedure: state.proc.id, solve_for: state.solve, options: { power_curve: false }, inputs };
  const form = state.form;
  clearResult();
  setBusy(true, state.proc.simulation ? 'Simulating…' : 'Calculating…');
  // the answer goes only under the form it was computed from (not one opened or reset since)
  const asked = () => state.form === form && state.solve === request.solve_for;
  for (let attempt = 0; ; attempt++) {
    try {
      const { result, engineMs } = await engine.run(request, MAX_ROWS);
      setBusy(false, '');
      if (asked()) renderResult(result, engineMs);
    } catch (e) {
      // the multi-core engine stopped: once more, on the one-thread engine that replaced it
      if (e.crashed && e.retry && attempt === 0 && asked()) {
        setBusy(true, 'Calculating again with one thread…');
        continue;
      }
      setBusy(false, e.stopped ? 'Stopped.' : '');
      if (!e.stopped && asked()) showError(e.message);
    }
    return;
  }
}

// ---------------------------------------------------------------------------------------------------------
// Speed test: every method and solve-for choice with its default inputs (simulation methods by their quick
// approximation), then three full simulations of 1,000 studies
// ---------------------------------------------------------------------------------------------------------
const FULL_SIMULATIONS = ['sim.means.robust', 'sim.lmm.cluster', 'sim.winratio'];

function percentile(sorted, q) {
  return sorted[Math.min(sorted.length - 1, Math.floor(q * sorted.length))];
}

async function speedTest() {
  if (state.busy || !state.catalog) return;
  const out = $('bench-out');
  setBusy(true);
  $('bench-stop').hidden = false;
  try {
    const requests = [];
    for (const p of state.catalog.procedures) {
      for (const s of (await schema(p.id)).solve_for) {
        requests.push({ name: p.name, request: { procedure: p.id, solve_for: s.value, options: { power_curve: true, approximate: true }, inputs: {} } });
      }
    }
    const times = [];
    const started = performance.now();
    for (const [i, r] of requests.entries()) {
      out.textContent = `Calculation ${i + 1} of ${requests.length}…`;
      const { engineMs } = await engine.run(r.request);
      times.push({ name: r.name, t: engineMs });
    }
    const wall = performance.now() - started;
    const sorted = times.map((x) => x.t).sort((a, b) => a - b);
    const slowest = times.reduce((a, b) => (b.t > a.t ? b : a));
    const sims = [];
    for (const id of FULL_SIMULATIONS) {
      const name = state.catalog.procedures.find((p) => p.id === id)?.name ?? id;
      out.textContent = `Simulating 1,000 studies: ${name}…`;
      const { result, engineMs } = await engine.run({ procedure: id, solve_for: 'power', options: { power_curve: false }, inputs: { reps: '1000' } });
      sims.push({ name, t: engineMs, threads: result.threads, ok: result.ok });
    }
    out.replaceChildren(
      el('table', { class: 'bench' },
        el('thead', {}, el('tr', {}, el('th', { scope: 'col' }, 'What'), el('th', { scope: 'col' }, 'Time'))),
        el('tbody', {},
          el('tr', {}, el('td', {}, `${requests.length} calculations: every method and solve-for choice, with a power curve`), el('td', {}, `${ms(wall)} in all`)),
          el('tr', {}, el('td', { class: 'sub' }, 'Median calculation'), el('td', {}, ms(percentile(sorted, 0.5)))),
          el('tr', {}, el('td', { class: 'sub' }, '9 in 10 take at most'), el('td', {}, ms(percentile(sorted, 0.9)))),
          el('tr', {}, el('td', { class: 'sub' }, `Slowest: ${slowest.name}`), el('td', {}, ms(slowest.t))),
          sims.map((s) => el('tr', {},
            el('td', {}, `Full simulation, 1,000 studies: ${s.name}${s.threads > 1 ? ` (${s.threads} threads)` : ''}`),
            el('td', {}, s.ok ? ms(s.t) : 'failed'))))));
  } catch (e) {
    out.textContent = e.stopped ? 'Stopped.' : `The speed test failed: ${e.message}`;
  } finally {
    $('bench-stop').hidden = true;
    setBusy(false);
  }
}

// ---------------------------------------------------------------------------------------------------------
async function boot() {
  $('form').addEventListener('submit', calculate);
  $('stop').addEventListener('click', () => engine.stop());
  $('activity-stop').addEventListener('click', () => engine.stop());
  $('bench-stop').addEventListener('click', () => engine.stop());
  $('bench').addEventListener('click', speedTest);
  $('search').addEventListener('input', onSearch);
  $('skip').addEventListener('click', () => {
    const target = !$('method').hidden ? $('m-name') : mode === 'guide' ? $('guide-questions').querySelector('input') : $('search');
    (target ?? $('mode-guide')).focus();
  });
  $('mode-guide').addEventListener('click', () => { openSeq++; setMode('guide'); });
  $('mode-browse').addEventListener('click', () => {
    openSeq++;
    setMode('browse');
    if (!state.proc && state.catalog) choose('means.two.t');
  });
  const edited = () => {
    state.form++;
    clearResult();
    if (guidedForm && !guide?.matches(state.proc, value)) {
      guide?.clearApplied();
      $('m-source').textContent = 'Calculator settings changed. Reapply the recommended analysis above to restore its choices.';
    }
  };
  $('fields').addEventListener('input', edited);
  $('fields').addEventListener('change', edited);
  engine.onstart = () => showStatus();
  window.addEventListener('online', () => engine.retry());
  window.addEventListener('hashchange', () => {
    const selection = parseSelectionHash(location.hash);
    if (known(selection.procedure)) {
      navigate(selection, false, false);
    } else if (location.hash.length > 1 && state.proc && !$('method').hidden) {
      // a malformed or unknown name: the method shown stays, and so does its address
      address(`${location.pathname}${location.search}#${state.proc.id}`);
    } else {
      openSeq++;
      state.form++;
      clearResult();
      guidedForm = false;
      setMode('guide');
    }
  });
  if (!(await showStatus())) {
    $('guide-questions').textContent = 'The study guide needs the calculator engine. Reload the page when the connection is available.';
    return;
  }
  try {
    state.catalog = (await engine.call('catalog')).value;
    for (const [id, proc] of Object.entries((await engine.call('procedures')).value)) schemas.set(id, proc);
    if ($('search').value.trim()) search();
    else renderList(null);
    guide = mountGuide($('guide'), {
      schemas,
      onChange() {
        openSeq++;
        guidedForm = false;
        state.form++;
        clearResult();
        $('method').hidden = true;
        $('navigation-error').hidden = true;
        address(location.pathname + location.search);
      },
      onApply: (selection) => navigate(selection, true),
    });
    const first = parseSelectionHash(location.hash);
    if (known(first.procedure)) await navigate(first, false, false);
    else if (mode === 'browse') {
      // the first method, in place of the bare address (not a step of its own in the browser's history)
      const start = { procedure: 'means.two.t', inputs: {} };
      if (await openMethod(start)) {  // not navigate(): no focus, so the page stays at its top
        setMode('browse');
        address(location.pathname + location.search + selectionHash(start));
      }
    }
    $('bench').setAttribute('aria-disabled', 'false');
  } catch (e) {
    $('engine-status').textContent = `The engine did not load: ${e.message}`;
  }
}

// The service worker keeps the page working offline and, on hosts that cannot send the headers, lets the engine
// use every core from the next visit on.  Where it is not allowed, the page works without it.
if ('serviceWorker' in navigator) navigator.serviceWorker.register('sw.js').catch(() => {});

boot();
