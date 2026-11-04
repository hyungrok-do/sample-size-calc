(() => {
// The genuine Lean parser, elaborator, interpreter and kernel run in this worker.
// Runtime initialization follows cauli/lean4-wasm-in-browser (Apache-2.0).
let manifest, cache, output = [], diagnostics = [], outputBytes = 0;
const base = new URL('lean-wasm/', self.location.href);
const progress = text => self.postMessage({ type: 'progress', text });
const hash = async bytes => [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))].map(x => x.toString(16).padStart(2, '0')).join('');
async function asset(file, expected) {
  const url = new URL(file, base); url.searchParams.set('v', manifest.version);
  let response = await cache?.match(url.href);
  if (!response) {
    response = await fetch(url);
    if (!response.ok) throw Error(`Could not load compiler asset ${file} (${response.status}).`);
    const bytes = new Uint8Array(await response.arrayBuffer());
    if (await hash(bytes) !== expected) throw Error(`Compiler asset checksum failed: ${file}`);
    // Keep the MIME type: cached JavaScript must still load as a module offline.
    try { await cache?.put(url.href, new Response(bytes, { headers: { 'Content-Type': response.headers.get('Content-Type') || 'application/octet-stream' } })); } catch {}
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
  if (/^\s*\[(?:DEBUG|WASM DEBUG)/.test(text)) {
    const match = /\[DEBUG:PROGRESS\]\s*(\d+)\/(\d+)/.exec(text);
    if (match) progress(`Loading Lean modules: ${match[1]} / ${match[2]}`);
    return;
  }
  outputBytes += text.length;
  if (outputBytes > 8 * 1024 * 1024) return;
  if (text.startsWith('{')) {
    try { const item = JSON.parse(text); if (item.severity) { diagnostics.push(item); return; } } catch {}
  }
  output.push(text);
}
function leanString(text) {
  const ptr = Module.stringToNewUTF8(text), result = Module._lean_mk_string(ptr);
  Module._free(ptr); return result;
}
function run(source) {
  output = []; diagnostics = []; outputBytes = 0;
  const result = Module._lean_wasm_compile(leanString(source), leanString('/workspace/source.lean'));
  if ((Module.getValue(result + 7, 'i8') & 255) !== 0) {
    Module._lean_io_result_show_error(result);
    throw Error(output.join('\n') || 'Lean could not import the requested module.');
  }
  const errors = diagnostics.filter(d => d.severity === 'error' || d.kind === 'hasSorry' || /sorryAx|declaration uses [`'"]sorry/.test(d.data));
  if (errors.length) throw Error(errors.map(d => `source.lean:${d.pos?.line || 1}:${(d.pos?.column || 0) + 1}: ${d.data}`).join('\n'));
  if (outputBytes > 8 * 1024 * 1024) throw Error('Compiler output exceeded 8 MB.');
  return { output, diagnostics };
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
  cache = await caches.open(`ssc-lean-${manifest.version}`).catch(() => null);
  progress(`Loading the browser compiler (${Math.ceil(manifest.downloadBytes / 1048576)} MB on first use)…`);
  const binary = await inflate(await asset('lean.wasm.pack', manifest.runtime['lean.wasm.pack'].sha256));
  // Authenticate the exact glue too. importScripts uses the same immutable URL.
  await asset('lean.js', manifest.runtime['lean.js'].sha256);
  const script = new URL('lean.js', base); script.searchParams.set('v', manifest.version);
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
  run('import SSC\n#check SSC.Report.exportReport');
  self.postMessage({ type: 'ready', version: manifest.leanVersion });
}
self.onmessage = async ({ data }) => {
  try {
    if (data.type === 'initialize') { await initialize(); return; }
    if (data.type !== 'compile') return;
    const source = data.source;
    if (typeof source !== 'string' || !source.trim() || new TextEncoder().encode(source).length > 131072) throw Error('Supply 1–131072 bytes of Lean source.');
    // The editor contains a complete ordinary Lean module, including its imports.
    // run_cmd is standard Lean syntax. A per-run token separates the authenticated
    // reporter's output from ordinary #eval output in the same stdout stream.
    const nonce = crypto.randomUUID();
    const prefix = `SSC_PROOF_REPORT ${nonce} `;
    const result = run(source + `\n\nrun_cmd _root_.SSC.Report.exportReport "${nonce}"\n`);
    // run_cmd uses Lean's output-capture machinery. The WASM frontend can
    // deliver its output as an information diagnostic instead of raw stdout.
    const lines = [...result.output, ...result.diagnostics.filter(d => d.severity === 'information').map(d => d.data)]
      .flatMap(text => text.split('\n'));
    const internal = line => line.startsWith(prefix) || line.startsWith('SSC_PROOF_STEP ') || line.startsWith('SSC_PROOF_TRACE_LIMIT ');
    const reports = lines.filter(text => text.startsWith(prefix));
    if (reports.length !== 1) throw Error('No checked declaration report. Include import SSC at the start of the Lean module.');
    const report = JSON.parse(reports[0].slice(prefix.length));
    const traces = lines.filter(text => text.startsWith('SSC_PROOF_STEP '))
      .map(text => JSON.parse(text.slice('SSC_PROOF_STEP '.length)));
    const declarations = report.declarations.sort((a, b) => a.line - b.line || a.column - b.column);
    for (const declaration of declarations) {
      const seen = new Set();
      declaration.traceTruncated = lines.includes('SSC_PROOF_TRACE_LIMIT ' + declaration.name);
      declaration.steps = traces.filter(step => step.declaration === declaration.name && !/^by(?:\s|$)/.test(step.tactic))
        .filter(step => { const key = JSON.stringify(step); if (seen.has(key)) return false; seen.add(key); return true; })
        .sort((a, b) => a.line - b.line || a.column - b.column);
    }
    self.postMessage({ type: 'result', id: data.id, value: {
      ok: true, compiler: 'lean-wasm-proofs-v2', language: 'Lean 4', leanVersion: manifest.leanVersion,
      packageVersion: manifest.version, reporterSourceSha256: manifest.sourceSha256,
      traceTruncated: lines.some(line => line.startsWith('SSC_PROOF_TRACE_LIMIT ')),
      results: declarations, module: source,
      output: result.output.flatMap(text => text.split('\n')).filter(line => !internal(line)).join('\n'),
      diagnostics: result.diagnostics.map(d => ({ ...d, data: d.data.split('\n').filter(line => !internal(line)).join('\n') }))
        .filter(d => d.data.trim()),
    } });
  } catch (error) {
    self.postMessage({ type: 'error', id: data.id, message: error?.message || 'Lean could not complete this compilation. Reload the compiler to retry.' });
  }
};
})();
