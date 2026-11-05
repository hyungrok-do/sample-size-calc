(() => {
// The genuine Lean parser, elaborator, interpreter and kernel run in this worker.
// Runtime initialization follows cauli/lean4-wasm-in-browser (Apache-2.0).
let manifest, cache, checker, output = [], diagnostics = [], outputBytes = 0;
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
  else if (count < 0) Atomics.sub(Module.HEAP32, value >>> 2, 1);
  return value;
}
function ioValue(result) {
  if ((Module.getValue(result + 7, 'i8') & 255) !== 0) {
    Module._lean_io_result_show_error(result);
    throw Error(output.join('\n') || 'Lean could not process the requested module.');
  }
  return Module.getValue(result + 8, 'i32');
}
function run(source, nonce) {
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
    const session = Module._lean_apply_3(checker, leanString(source), leanString(nonce), 1);
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
  const drivers = [...bootstrap.output, ...bootstrap.diagnostics.filter(d => d.severity === 'information').map(d => d.data)]
    .flatMap(text => text.split('\n')).filter(text => /^SSC_BROWSER_DRIVER \d+$/.test(text));
  if (drivers.length !== 1) throw Error('Missing trusted compiler runtime entry point.');
  checker = Number(drivers[0].slice('SSC_BROWSER_DRIVER '.length));
  // LeanClosure = 245, with reference count 0 for Runtime.markPersistent.
  if (!Number.isSafeInteger(checker) || checker <= 0 || checker % 8 !== 0
      || (Module.getValue(checker + 7, 'i8') & 255) !== 245 || Module.getValue(checker, 'i32') !== 0)
    throw Error('Incompatible compiler runtime entry point.');
  // A ready compiler must also be ready for the immediately-following offline
  // reload. Writes overlap validation/decompression/imports, but finish here.
  await Promise.all(cacheWrites);
  self.postMessage({ type: 'ready', version: manifest.leanVersion });
}
self.onmessage = async ({ data }) => {
  try {
    if (data.type === 'initialize') { await initialize(); return; }
    if (data.type !== 'compile') return;
    const compileStarted = performance.now();
    const source = data.source;
    if (typeof source !== 'string' || !source.trim() || new TextEncoder().encode(source).length > 131072) throw Error('Supply 1–131072 bytes of Lean source.');
    // User syntax is parsed by Lean inside a trusted driver. Neither the nonce
    // nor the final audit command is appended to the user-controlled module.
    const nonce = crypto.randomUUID();
    const prefix = 'SSC_COMPILE_RESULT ' + nonce + ' ';
    const result = run(source, nonce);
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
    self.postMessage({ type: 'result', id: data.id, value: {
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
    const fatal = !(error instanceof Error) || error instanceof RangeError || error instanceof WebAssembly.RuntimeError;
    self.postMessage({ type: 'error', id: data.id, fatal, message: error?.message || 'Lean execution stopped unexpectedly. The compiler will restart for the next request.' });
  }
};
})();
