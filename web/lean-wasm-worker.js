(() => {
// The genuine Lean parser, elaborator, interpreter and kernel run in this worker.
// Runtime initialization follows cauli/lean4-wasm-in-browser (Apache-2.0).
let manifest, cache, numeric, cppRuntime, output = [], diagnostics = [], outputBytes = 0;
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
  if (errors.length) throw Error(errors.map(d => `source.lean:${Math.max(1, (d.pos?.line || 1) - 1)}: ${d.data}`).join('\n'));
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
  if (manifest.schema !== 1 || !Array.isArray(manifest.packs)) throw Error('Unsupported browser compiler package.');
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
  progress('Loading the C++ WebAssembly numerical backend…');
  // All four files belong to the same checksum-pinned static package.
  const numericalBytes = await asset('study.wasm', manifest.runtime['study.wasm'].sha256);
  await asset('study.mjs', manifest.runtime['study.mjs'].sha256);
  const numericalURL = new URL('study.mjs', base); numericalURL.searchParams.set('v', manifest.version);
  const factory = (await import(numericalURL.href)).default;
  const numericalModule = await WebAssembly.compile(numericalBytes);
  numeric = await factory({ instantiateWasm(imports, receive) {
    const instance = new WebAssembly.Instance(numericalModule, imports);
    receive(instance, numericalModule); return instance.exports;
  } });
  cppRuntime = new TextDecoder().decode(await asset('study-runtime.txt', manifest.runtime['study-runtime.txt'].sha256));
  progress('Initializing the statistical Lean compiler…');
  run('import SSC\n#check DesignCompiler.Planning.Study');
  self.postMessage({ type: 'ready', version: manifest.leanVersion });
}
self.onmessage = async ({ data }) => {
  try {
    if (data.type === 'initialize') { await initialize(); return; }
    if (data.type !== 'compile') return;
    const source = data.source;
    if (typeof source !== 'string' || !source.trim() || new TextEncoder().encode(source).length > 131072) throw Error('Supply 1–131072 bytes of Lean source.');
    const module = 'import SSC\n' + source + '\n';
    const result = run(module + '\n#export_compiler');
    const reports = result.output.filter(text => text.startsWith('SSC_LEAN_REPORT '));
    if (reports.length !== 1) throw Error('Lean did not produce a unique checked study report.');
    const plans = JSON.parse(reports[0].slice('SSC_LEAN_REPORT '.length));
    if (!Array.isArray(plans) || !plans.length) throw Error('Add a solve sample_size, power, or effect command.');
    for (const plan of plans) { const lowered = lower(plan); plan.evaluation = evaluate(lowered); plan.cpp = emitCpp(lowered); }
    self.postMessage({ type: 'result', id: data.id, value: { ok: true, compiler: 'lean-wasm-statistical-v1', language: 'Lean 4',
      leanKernelChecked: true, results: plans, module,
      diagnostics: result.diagnostics.map(d => d.data).join('\n') } });
  } catch (error) {
    self.postMessage({ type: 'error', id: data.id, message: error?.message || 'Lean could not complete this compilation. Reload the compiler to retry.' });
  }
};

function lower(plan) {
  if (plan.kind !== 'statistical_design') throw Error('The browser compiler expects a statistical study.');
  const rat = value => {
    const [a,b='1'] = String(value).split('/'), result = Number(a)/Number(b);
    if (!Number.isFinite(result)) throw Error('Planning values exceed numerical precision. Rescale the study.');
    return result;
  };
  const count = value => {
    if (!Number.isSafeInteger(value) || value < 0 || value > 10000000) throw Error('Invalid planning count.');
    return value;
  };
  return { variance:rat(plan.statistics.variance), effect:rat(plan.statistics.effect),
    alpha:rat(plan.criterion.alpha), target:rat(plan.criterion.targetPower),
    minimum:count(plan.criterion.minN), maximum:count(plan.criterion.maxN), fixed:count(plan.fixedN), mode:plan.solveFor,
    groups:plan.groups.map(g=>({label:g.label, units:count(g.unitsPerN), participants:count(g.participantsPerUnit)})) };
}
function evaluate(plan) {
  const result = JSON.parse(numeric.cwrap('ssc_evaluate_study', 'string', ['string'])(JSON.stringify(plan)));
  if (result.error) throw Error(result.error);
  return result;
}
function cppString(value) {
  const bytes = new TextEncoder().encode(value);
  const text = [...bytes].map(c => c >= 32 && c < 127 && c !== 34 && c !== 92 ? String.fromCharCode(c) : '\\' + c.toString(8).padStart(3,'0')).join('');
  return `std::string("${text}", ${bytes.length})`;
}
function emitCpp(p) {
  const groups = p.groups.map(g=>`{${cppString(g.label)}, ${g.units}, ${g.participants}}`).join(',');
  return '// Generated from a Lean-checked statistical study. The same numerical runtime executes in browser WebAssembly.\n'
    + cppRuntime.replace('#pragma once','') + '\n#include <iostream>\nint main() {\n  try {\n'
    + `    const ssc_lean::Plan plan{${p.variance.toPrecision(17)}, ${p.effect.toPrecision(17)}, ${p.alpha.toPrecision(17)}, ${p.target.toPrecision(17)}, ${p.minimum}, ${p.maximum}, ${p.fixed}, ${cppString(p.mode)}, {${groups}}};\n`
    + '    std::cout << ssc_lean::evaluate(plan) << "\\n";\n    return 0;\n'
    + '  } catch(const std::exception &e) { std::cerr << e.what() << "\\n"; return 2; }\n}\n';
}

})();
