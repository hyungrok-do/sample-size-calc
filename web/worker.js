// The calculator's engine (the C++ core compiled to WebAssembly) in a Web Worker, so a long calculation never
// freezes the page.  Messages: {id, op, arg} in, {id, ok, value, ms} or {id, ok: false, error} out, and one
// {stage: 'compiled'} once the engine is downloaded and compiled (what is left of starting is quick).
const started = performance.now();
const variant = new URL(self.location.href).searchParams.get('engine') === 'ssc_web_mt' ? 'ssc_web_mt' : 'ssc_web';
let loadMs = 0;

async function compile() {
  const response = await fetch(new URL(`./engine/${variant}.wasm?v=dev`, import.meta.url));
  if (!response.ok) throw new Error(`the engine file could not be downloaded (${response.status})`);
  const type = response.headers.get('content-type') ?? '';
  // a static host that answers a missing file with its home page
  if (type.includes('text/html')) throw new Error('the engine file is missing (the server sent a web page instead)');
  // a server that sends another type for .wasm: compile from the bytes
  return type.includes('application/wasm')
    ? WebAssembly.compileStreaming(response)
    : WebAssembly.compile(await response.arrayBuffer());
}

async function load() {
  if (variant === 'ssc_web_mt' && !(self.crossOriginIsolated && typeof SharedArrayBuffer === 'function')) {
    throw new Error('the multi-core engine needs a cross-origin isolated page');
  }
  const [module, compiled] = await Promise.all([import(`./engine/${variant}.mjs?v=dev`), compile()]);
  self.postMessage({ stage: 'compiled' });
  const M = await module.default({
    instantiateWasm(imports, receive) {
      WebAssembly.instantiate(compiled, imports).then((instance) => receive(instance, compiled));
      return {};
    },
  });
  loadMs = performance.now() - started;
  return M;
}

const engine = load();
engine.catch(() => {});  // reported to each call below

// A JavaScript string in the engine's memory, for the length of `use`.
function withText(M, text, use) {
  const size = M.lengthBytesUTF8(text) + 1;
  const ptr = M._malloc(size);
  M.stringToUTF8(text, ptr, size);
  try {
    return use(ptr);
  } finally {
    M._free(ptr);
  }
}

// A string the engine returned, freed after reading.
function take(M, ptr) {
  if (!ptr) return null;
  const text = M.UTF8ToString(ptr);
  M._ssc_free(ptr);
  return text;
}

// How the engine's .wasm file arrived: bytes over the network (few or none from the cache), compressed and as run.
function download() {
  const entry = performance.getEntriesByType('resource').find((e) => e.name.includes('.wasm'));
  return entry
    ? { transfer: entry.transferSize, encoded: entry.encodedBodySize, decoded: entry.decodedBodySize, ms: entry.duration }
    : null;
}

const procedure = (M, id) => withText(M, id, (p) => JSON.parse(take(M, M._ssc_procedure_json(p))));

const ops = {
  info: (M) => ({
    version: M.UTF8ToString(M._ssc_version()),
    threads: M._ssc_hardware_threads(),
    variant,
    loadMs,
    download: download(),
  }),
  catalog: (M) => JSON.parse(take(M, M._ssc_catalog_json())),
  procedure,
  // every procedure's description at once, so that opening a method never waits for a calculation
  procedures: (M) => {
    const all = {};
    for (const p of JSON.parse(take(M, M._ssc_catalog_json())).procedures) all[p.id] = procedure(M, p.id);
    return all;
  },
  search: (M, { query, limit }) => withText(M, query, (p) => JSON.parse(take(M, M._ssc_search_json(p, limit)))),
  // maxRows: the rows a page shows (up to 20,000 scenarios would take seconds to copy to it and draw)
  run: (M, { request, maxRows = Infinity }) => {
    const token = M._ssc_token_new();
    try {
      return withText(M, JSON.stringify(request), (p) => {
        const t = performance.now();
        const text = take(M, M._ssc_run_json(p, token));
        const engineMs = performance.now() - t;
        const result = JSON.parse(text);
        delete result.plots;  // the page draws no charts
        if (result.rows?.length > maxRows) {
          result.rows_total = result.rows.length;
          result.rows = result.rows.slice(0, maxRows);
          if (result.dropout?.rows) result.dropout.rows = result.dropout.rows.filter((r) => r._row < maxRows);
        }
        // the engine's memory never shrinks: the page replaces an engine that a large answer left big
        return { result, engineMs, memoryBytes: M.HEAP8.buffer.byteLength };
      });
    } finally {
      M._ssc_token_free(token);
    }
  },
};

self.onmessage = async (event) => {
  const { id, op, arg } = event.data ?? {};
  if (id === undefined) return;
  try {
    const M = await engine;
    if (!ops[op]) throw new Error(`unknown request '${op}'`);
    const t = performance.now();
    const value = ops[op](M, arg);
    self.postMessage({ id, ok: true, value, ms: performance.now() - t });
  } catch (e) {
    self.postMessage({ id, ok: false, error: String(e?.message ?? e) });
  }
};
