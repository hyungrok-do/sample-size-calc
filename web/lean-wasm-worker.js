(() => {
// The genuine Lean parser, elaborator, interpreter and kernel run in this worker.
// Runtime initialization follows cauli/lean4-wasm-in-browser (Apache-2.0).
let manifest, cache, checker, evaluator, prepared = null, output = [], diagnostics = [], outputBytes = 0;
const cacheWrites = new Set();
const base = new URL('lean-wasm/', self.location.href);
const progress = text => self.postMessage({ type: 'progress', text });
const hash = async bytes => [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))].map(x => x.toString(16).padStart(2, '0')).join('');
async function asset(file, expected) {
  // Unchanged runtime/library content survives a new site or reporter build.
  // Authentication still runs on every cache read before any bytes are used.
  const url = new URL(file, base); url.searchParams.set('v', expected);
  let response = await cache?.match(url.href);
  if (!response) {
    response = await fetch(url);
    if (!response.ok) throw Error(`Could not load compiler asset ${file} (${response.status}).`);
    const bytes = new Uint8Array(await response.arrayBuffer());
    if (await hash(bytes) !== expected) throw Error(`Compiler asset checksum failed: ${file}`);
    // Keep the MIME type: cached JavaScript must still load as a module offline.
    if (cache) {
      const write = cache.put(url.href, new Response(bytes, { headers: { 'Content-Type': response.headers.get('Content-Type') || 'application/octet-stream' } }))
        .catch(() => {}).finally(() => cacheWrites.delete(write));
      cacheWrites.add(write);
    }
    return bytes;
  }
  const bytes = new Uint8Array(await response.arrayBuffer());
  if (await hash(bytes) !== expected) { await cache.delete(url.href); throw Error(`Cached compiler asset is damaged: ${file}. Reload the compiler.`); }
  return bytes;
}
async function inflate(bytes) {
  return new Uint8Array(await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());
}
function mkdir(path) {
  let current = '';
  for (const part of path.split('/').filter(Boolean)) {
    current += '/' + part;
    if (!Module.FS.analyzePath(current).exists) Module.FS.mkdir(current);
  }
}
function line(text) {
  outputBytes += new TextEncoder().encode(text).byteLength;
  if (outputBytes > 8 * 1024 * 1024) return;
  if (/^\s*\[(?:DEBUG|WASM DEBUG)/.test(text)) {
    const match = /\[DEBUG:PROGRESS\]\s*(\d+)\/(\d+)/.exec(text);
    if (match) progress(`Loading Lean modules: ${match[1]} / ${match[2]}`);
    return;
  }
  if (text.startsWith('{')) {
    try { const item = JSON.parse(text); if (item.severity) { diagnostics.push(item); return; } } catch {}
  }
  output.push(text);
}
function leanString(text) {
  // Use the length-aware constructor: a NUL inside a source comment must not
  // truncate the rest of the module at the C-string boundary.
  const ptr = Module.stringToNewUTF8(text);
  const result = Module._lean_mk_string_from_bytes(ptr, new TextEncoder().encode(text).byteLength);
  Module._free(ptr); return result;
}
function release(result) {
  // Pinned Lean runtime ABI: mirror lean_dec_ref from lean/lean.h. The cold
  // runtime function handles the final reference and multi-threaded objects.
  const count = Module.getValue(result, 'i32');
  if (count > 1) Module.setValue(result, count - 1, 'i32');
  else if (count !== 0) Module._lean_dec_ref_cold(result);
}
function retain(value) {
  const count = Module.getValue(value, 'i32');
  if (count > 0) Module.setValue(value, count + 1, 'i32');
  else if (count < 0) Atomics.sub(new Int32Array(Module.wasmMemory.buffer), value >>> 2, 1);
  return value;
}
function ioValue(result) {
  if ((Module.getValue(result + 7, 'i8') & 255) !== 0) {
    Module._lean_io_result_show_error(result);
    throw Error(output.join('\n') || 'Lean could not process the requested module.');
  }
  return Module.getValue(result + 8, 'i32');
}
function readLeanString(value) {
  // Pinned wasm32 Lean ABI: header (8), size/capacity/UTF8 length (4 each),
  // then bytes. Size includes the terminal NUL; decode by length, not sentinel.
  if ((Module.getValue(value + 7, 'i8') & 255) !== 249) throw Error('Invalid Lean string result.');
  const size = Module.getValue(value + 8, 'i32') >>> 0;
  if (size < 1 || size - 1 + outputBytes > 8 * 1024 * 1024) throw Error('Compiler output exceeded 8 MB.');
  // This pinned glue exports the supplied memory, not Emscripten heap views.
  // Recreate the view after each call because the runtime can grow memory.
  const heap = new Uint8Array(Module.wasmMemory.buffer);
  if (value + 20 + size > heap.length) throw Error('Invalid Lean string size.');
  return new TextDecoder('utf-8', { fatal: true }).decode(heap.slice(value + 20, value + 20 + size - 1));
}
function validateJsonNumber(value) {
  if (typeof value === 'number' && (!Number.isFinite(value)
      || (Number.isInteger(value) && !Number.isSafeInteger(value))))
    throw Error('JSON numbers must be finite and integers must be within the safe integer range. Encode exact integers and rationals as strings.');
  return value;
}
function runEvaluation(input) {
  output = []; diagnostics = []; outputBytes = 0;
  const started = performance.now();
  const result = Module._lean_apply_2(evaluator, leanString(input), 1);
  try {
    const text = readLeanString(ioValue(result));
    const envelope = JSON.parse(text, (_key, value) => validateJsonNumber(value));
    if (!envelope.ok) throw Error(envelope.error || 'Prepared evaluation failed.');
    return { value: envelope.value, evaluationMs: performance.now() - started };
  } finally { release(result); output = []; diagnostics = []; outputBytes = 0; }
}
function run(source, nonce, entry = '') {
  const started = performance.now();
  const timing = { prepareMs: 0, commandMs: 0, driverMs: 0 };
  output = []; diagnostics = []; outputBytes = 0;
  if (nonce === undefined) {
    const result = Module._lean_wasm_compile(leanString(source), leanString('/workspace/source.lean'));
    try { ioValue(result); } finally { release(result); }
  } else {
    // Each returned action is Lean's native frontend specialized to the current
    // command. Execute it with a flat host stack, outside the interpreted driver.
    // IO results own their fields; retain before passing a consumed argument.
    const session = Module._lean_apply_4(checker, leanString(source), leanString(nonce), leanString(entry), 1);
    timing.prepareMs = performance.now() - started;
    try {
      const next = ioValue(session);
      let count = 0, lastProgress = started;
      for (;;) {
        const driverStarted = performance.now();
        const step = Module._lean_apply_1(retain(next), 1);
        timing.driverMs += performance.now() - driverStarted;
        try {
          const option = ioValue(step);
          if (option === 1) break; // Option.none
          const action = Module.getValue(option + 8, 'i32');
          const commandStarted = performance.now();
          if (commandStarted - lastProgress >= 500) {
            progress(`Processing Lean command ${count + 1} (${((commandStarted - started) / 1000).toFixed(1)} s)…`);
            lastProgress = commandStarted;
          }
          const result = Module._lean_apply_1(retain(action), 1);
          try { ioValue(result); } finally { release(result); }
          timing.commandMs += performance.now() - commandStarted;
          count++;
        } finally { release(step); }
        if (outputBytes > 8 * 1024 * 1024) throw Error('Compiler output exceeded 8 MB.');
      }
    } finally { release(session); }
  }
  const errors = diagnostics.filter(d => d.severity === 'error' || d.kind === 'hasSorry'
    || (d.severity === 'warning' && /declaration uses [`'"]sorry/.test(d.data)));
  if (errors.length) throw Error(errors.map(d => `source.lean:${d.pos?.line || 1}:${(d.pos?.column || 0) + 1}: ${d.data}`).join('\n'));
  if (outputBytes > 8 * 1024 * 1024) throw Error('Compiler output exceeded 8 MB.');
  return { output, diagnostics, timing: { ...timing, totalMs: performance.now() - started } };
}
async function initialize() {
  let response;
  const metadataCache = await caches.open('ssc-lean-metadata').catch(() => null);
  const manifestURL = new URL('manifest.json', base);
  try {
    response = await fetch(manifestURL, { cache: 'no-cache' });
    if (!response.ok) throw Error('The browser compiler has not been packaged with this site.');
    try { await metadataCache?.put(manifestURL.href, response.clone()); } catch {}
  } catch (error) {
    response = await metadataCache?.match(manifestURL.href);
    if (!response) throw error;
  }
  manifest = await response.json();
  if (manifest.schema !== 2 || !Array.isArray(manifest.packs)) throw Error('Unsupported browser compiler package.');
  cache = await caches.open('ssc-lean-content-v1').catch(() => null);
  progress(`Loading the browser compiler (${Math.ceil(manifest.downloadBytes / 1048576)} MB on first use)…`);
  const binary = await inflate(await asset('lean.wasm.pack', manifest.runtime['lean.wasm.pack'].sha256));
  // Authenticate the exact glue too. importScripts uses the same immutable URL.
  await asset('lean.js', manifest.runtime['lean.js'].sha256);
  const script = new URL('lean.js', base); script.searchParams.set('v', manifest.runtime['lean.js'].sha256);
  await new Promise((resolve, reject) => {
    self.Module = {
      wasmBinary: binary, noInitialRun: true, thisProgram: '/bin/lean',
      wasmMemory: new WebAssembly.Memory({ initial: 8192, maximum: 65536, shared: true }),
      INITIAL_MEMORY: 512 * 1024 * 1024,
      mainScriptUrlOrBlob: script.href,
      locateFile: file => new URL(file, base).href,
      print: line, printErr: line,
      preRun: [() => { for (const dir of ['/lib/lean', '/bin', '/workspace']) mkdir(dir); Module.ENV.LEAN_PATH = '/lib/lean'; Module.FS.chdir('/workspace'); }],
      onRuntimeInitialized: () => {
        try {
          Module._lean_initialize_runtime_module(); Module._lean_initialize(); Module._lean_io_mark_end_initialization();
          Module._lean_init_task_manager?.(); Module._lean_enable_initializer_execution?.();
          const result = Module._lean_init_search_path();
          if ((Module.getValue(result + 7, 'i8') & 255) !== 0) throw Error('Could not initialize Lean module paths.');
          Module._emscripten_runtime_keepalive_push?.(); resolve();
        } catch (error) { reject(error); }
      },
      onAbort: error => reject(Error(`Lean WebAssembly stopped: ${error}`)),
    };
    try { importScripts(script.href); } catch (error) { reject(error); }
  });
  let finished = 0;
  const queue = [...manifest.packs];
  await Promise.all(Array.from({ length: 3 }, async () => {
    while (queue.length) {
      const pack = queue.shift(), bytes = await inflate(await asset(pack.file, pack.sha256));
      if (bytes.length !== pack.bytes) throw Error(`Incorrect compiler pack size: ${pack.file}`);
      for (const entry of pack.entries) {
        if (!/^(?:Init|Std|Lean|SSC)(?:\.|\/)/.test(entry.path) || entry.path.split('/').includes('..') || entry.offset < 0 || entry.offset + entry.bytes > bytes.length)
          throw Error('Invalid compiler library entry.');
        const path = '/lib/lean/' + entry.path; mkdir(path.slice(0, path.lastIndexOf('/')));
        Module.FS.writeFile(path, bytes.subarray(entry.offset, entry.offset + entry.bytes));
      }
      progress(`Preparing compiler libraries: ${++finished} / ${manifest.packs.length}`);
    }
  }));
  progress('Initializing Lean proof checking…');
  const bootstrap = run('import SSC\n#ssc_browser_driver');
  const bootstrapLines = [...bootstrap.output, ...bootstrap.diagnostics.filter(d => d.severity === 'information').map(d => d.data)]
    .flatMap(text => text.split('\n'));
  function entryPoint(label) {
    const entries = bootstrapLines.filter(text => new RegExp('^' + label + ' \\d+$').test(text));
    if (entries.length !== 1) throw Error('Missing trusted compiler runtime entry point: ' + label);
    const pointer = Number(entries[0].slice(label.length + 1));
    // LeanClosure = 245, with reference count 0 for Runtime.markPersistent.
    if (!Number.isSafeInteger(pointer) || pointer <= 0 || pointer % 8 !== 0
        || (Module.getValue(pointer + 7, 'i8') & 255) !== 245 || Module.getValue(pointer, 'i32') !== 0)
      throw Error('Incompatible compiler runtime entry point: ' + label);
    return pointer;
  }
  checker = entryPoint('SSC_BROWSER_DRIVER');
  evaluator = entryPoint('SSC_BROWSER_EVALUATOR');
  // A ready compiler must also be ready for the immediately-following offline
  // reload. Writes overlap validation/decompression/imports, but finish here.
  await Promise.all(cacheWrites);
  self.postMessage({ type: 'ready', version: manifest.leanVersion });
}
async function handleRequest(data) {
  try {
    if (data.type === 'initialize') { await initialize(); return; }
    if (data.type === 'evaluate') {
      const started = performance.now();
      if (!prepared || data.preparedId !== prepared.id) throw Error('Prepared calculator is unavailable or stale. Prepare the source again.');
      const input = JSON.stringify(data.parameters, (_key, value) => {
        validateJsonNumber(value);
        if (value === undefined || typeof value === 'bigint') throw Error('Parameters must be JSON values.');
        return value;
      });
      if (typeof input !== 'string' || new TextEncoder().encode(input).length > 131072) throw Error('Supply at most 131072 bytes of JSON parameters.');
      const evaluation = runEvaluation(input);
      self.postMessage({ type: 'result', id: data.id, value: {
        ok: true, kind: 'prepared-evaluation', verification: 'compiled-execution',
        concreteResultKernelChecked: false, prepared, value: evaluation.value,
        performance: { evaluationMs: evaluation.evaluationMs, totalMs: performance.now() - started },
      } });
      return;
    }
    if (data.type !== 'compile' && data.type !== 'prepare') return;
    // Every source request invalidates the prior artifact, including malformed
    // input or failed preparation. Evaluation cannot reach a stale closure.
    prepared = null;
    runEvaluation(''); // Release the Lean closure even when input validation fails.
    const compileStarted = performance.now();
    const source = data.source;
    if (typeof source !== 'string' || !source.trim() || new TextEncoder().encode(source).length > 131072) throw Error('Supply 1–131072 bytes of Lean source.');
    // User syntax is parsed by Lean inside a trusted driver. Neither the nonce
    // nor the final audit command is appended to the user-controlled module.
    const nonce = crypto.randomUUID();
    const prefix = 'SSC_COMPILE_RESULT ' + nonce + ' ';
    const entry = data.type === 'prepare' ? data.entry : '';
    if (typeof entry !== 'string' || (data.type === 'prepare' && (!entry.trim() || new TextEncoder().encode(entry).length > 1024)))
      throw Error('Supply a Lean declaration name of 1–1024 bytes.');
    const result = run(source, nonce, entry);
    const lines = [...result.output, ...result.diagnostics.filter(d => d.severity === 'information').map(d => d.data)]
      .flatMap(text => text.split('\n'));
    const reports = lines.filter(text => text.startsWith(prefix));
    if (reports.length !== 1) throw Error('No independently checked declaration report.');
    const envelope = JSON.parse(reports[0].slice(prefix.length));
    const sourceDiagnostics = envelope.diagnostics || [];
    const errors = sourceDiagnostics.filter(d => d.severity === 'error' || d.kind === 'hasSorry'
      || (d.severity === 'warning' && /declaration uses [\u0060'"]sorry/.test(d.data)));
    if (!envelope.ok || errors.length) throw Error(envelope.error || errors.map(d =>
      'source.lean:' + (d.pos?.line || 1) + ':' + ((d.pos?.column || 0) + 1) + ': ' + d.data).join('\n')
      || 'Lean could not check this module.');
    const report = envelope.report;
    const traces = report.steps, limits = report.traceLimits;
    const declarations = report.declarations.sort((a, b) => a.line - b.line || a.column - b.column);
    // Index traces once. Scanning the whole trace for every declaration makes
    // report assembly quadratic in the size of a module.
    const byDeclaration = new Map();
    for (const step of traces) {
      if (/^by(?:\s|$)/.test(step.tactic)) continue;
      let entry = byDeclaration.get(step.declaration);
      if (!entry) byDeclaration.set(step.declaration, entry = { seen: new Set(), steps: [] });
      const key = JSON.stringify(step);
      if (!entry.seen.has(key)) { entry.seen.add(key); entry.steps.push(step); }
    }
    const truncated = new Set(limits);
    for (const declaration of declarations) {
      declaration.traceTruncated = truncated.has(declaration.name);
      declaration.steps = (byDeclaration.get(declaration.name)?.steps || [])
        .sort((a, b) => a.line - b.line || a.column - b.column);
    }
    if (data.type === 'prepare') prepared = {
      id: crypto.randomUUID(), sourceSha256: await hash(new TextEncoder().encode(source)),
      entry, packageVersion: manifest.version,
    };
    self.postMessage({ type: 'result', id: data.id, value: {
      ...(prepared ? { prepared } : {}),
      ok: true, compiler: 'lean-wasm-proofs-v2', language: 'Lean 4', leanVersion: manifest.leanVersion,
      packageVersion: manifest.version, reporterSourceSha256: manifest.sourceSha256,
      performance: { ...result.timing, ...envelope.performance, ...report.performance,
        totalMs: performance.now() - compileStarted },
      traceTruncated: limits.length > 0,
      results: declarations, module: source,
      output: lines.filter(line => !line.startsWith(prefix)).join('\n'),
      diagnostics: sourceDiagnostics,
    } });
  } catch (error) {
    // A host stack overflow or WASM trap can leave the runtime unusable.
    // Normal Lean diagnostics keep the worker; runtime failures require restart.
    let fatal = !(error instanceof Error) || error instanceof RangeError || error instanceof WebAssembly.RuntimeError;
    if (fatal || data.type === 'compile' || data.type === 'prepare') prepared = null;
    if (!fatal && evaluator && (data.type === 'compile' || data.type === 'prepare')) {
      try { runEvaluation(''); } catch { fatal = true; }
    }
    self.postMessage({ type: 'error', id: data.id, fatal, message: error?.message || 'Lean execution stopped unexpectedly. The compiler will restart for the next request.' });
  }
}
// Hashing and initialization yield to the event loop. Serialize requests so a
// preparation handle can never be paired with another source's closure.
let requests = Promise.resolve();
self.onmessage = ({ data }) => { requests = requests.then(() => handleRequest(data)); };
})();
