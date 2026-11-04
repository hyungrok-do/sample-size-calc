#!/usr/bin/env python3
"""Build the browser compiler using the actual pinned Lean WASM compiler.

No native Lean installation, server, or private Actions artifact is required.
Downloaded archives are immutable releases authenticated before extraction.
"""
import argparse
import ast
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
LEAN_COMMIT = '62b6a2291302d4bbeace37642a066b7510d0145c'
RELEASES = {
    'runtime': ('https://github.com/cauli/lean4-wasm-in-browser/releases/download/runtime-62b6a22-editor1/lean-runtime-fixture.tar.gz',
                'a66af27e859464db9849df10ea3bd7aa3b5ee824fbcc7c330fb8cf5931cc5527'),
    'static': ('https://github.com/cauli/lean4-wasm-in-browser/releases/download/runtime-62b6a22-compact1/pages-assets-62b6a22-compact1.tar.gz',
               'e0c1ae4d40fd5c4da18d4a38e254e8adef603b855757daff103659a4b80c9a39'),
}

def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()

def download(name, cache):
    url, sha = RELEASES[name]
    path = cache / (name + '.tar.gz')
    if path.exists() and digest(path) == sha:
        return path
    print(f'Downloading pinned {name} artifact', flush=True)
    temporary = path.with_suffix('.download')
    # curl handles GitHub release redirects/retries without any authentication.
    subprocess.run(['curl', '-fsSL', '--retry', '3', url, '-o', str(temporary)], check=True)
    if digest(temporary) != sha:
        raise RuntimeError(f'{name} artifact checksum mismatch')
    temporary.replace(path)
    return path

def unpack(archive, kit, runtime):
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            name = member.name.removeprefix('./')
            relative = None
            if runtime and name in ['bin/lean.js', 'bin/lean.wasm', 'bin/package.json']:
                relative = Path(name)
            elif not runtime and name.startswith('lean-lib/'):
                rel = name.removeprefix('lean-lib/')
                module = rel.split('/')[0].split('.')[0]
                if module in ['Init', 'Std', 'Lean'] and rel.endswith(('.olean', '.ir', '.ir.sig')):
                    relative = Path('lib/lean') / rel
            if relative is None:
                continue
            if not member.isfile() or '..' in relative.parts or relative.is_absolute():
                raise RuntimeError(f'Unexpected archive entry: {name}')
            dest = kit / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as source, dest.open('wb') as output:
                shutil.copyfileobj(source, output)

def browser_source():
    # The model and elaborator are the actual source used by the native compiler.
    # The browser imports Lean's rational arithmetic instead of the much larger
    # native Mathlib probability audit and separate general real-equation solver.
    model = (ROOT/'lean/DesignCompiler/Planning.lean').read_text()
    model = model[model.index('namespace DesignCompiler.Planning'):model.index('/-- The off-diagonal')]
    command = (ROOT/'lean/DesignCompiler/Planning/Command.lean').read_text()
    command = command[command.index('open Lean Meta Elab Command'):]
    substitutions = {
        '  let equations := (DesignCompiler.Equation.resultExtension.getState env).map toJson\n': '',
        '(studies ++ equations)': 'studies',
        '  let literal ← Term.elabTermAndSynthesize (← `(($(quote value.num) : ℚ) / ($(quote value.den) : ℚ))) (some (mkConst ``Rat))':
            '  let numerator ← if value.num < 0 then `(term| - ($(quote value.num.natAbs) : ℚ)) else `(term| ($(quote value.num.natAbs) : ℚ))\n'
            '  let literal ← Term.elabTermAndSynthesize (← `(($numerator) / ($(quote value.den) : ℚ))) (some (mkConst ``Rat))',
        'and the generic normal-tail identity.':
            'in this browser. The generic probability identities have a separate native Mathlib audit.',
        'Uses the checked shifted-normal rejection-probability identity.':
            'Uses the shifted-normal rejection-probability identity audited in the native Mathlib library.',
    }
    for before, after in substitutions.items():
        if command.count(before) != 1:
            raise RuntimeError(f'Browser compiler adapter needs review: {before}')
        command = command.replace(before, after)
    return ('import Lean\n\nnotation "ℚ" => Rat\n\nnamespace DesignCompiler\n'
            'def rationalJson (x : Rat) : Lean.Json := Lean.Json.str s!"{x.num}/{x.den}"\n'
            'end DesignCompiler\n\n' + model + '\nend DesignCompiler.Planning\n\n' + command)

def node_command():
    # emsdk_env.sh adds emcc to PATH, but can leave an older system Node first.
    # Use the SDK's configured engine, as emcc does. Older V8 versions reject
    # the genuine Lean runtime's 104769 exports before running any Lean code.
    command = ast.literal_eval(subprocess.check_output(['em-config', 'NODE_JS'], text=True).strip())
    if isinstance(command, str):
        command = [command]
    if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
        raise RuntimeError('Emscripten NODE_JS must name a JavaScript runtime')
    return command

def package(kit, out, source):
    out.mkdir(parents=True, exist_ok=True)
    # Compressed opaque extensions avoid a static host decompressing them twice.
    wasm = gzip.compress((kit/'bin/lean.wasm').read_bytes(), compresslevel=9, mtime=0)
    (out/'lean.wasm.pack').write_bytes(wasm)
    glue = (kit/'bin/lean.js').read_text()
    # This compact upstream glue omitted wasmBinary from its incoming API.
    # Restore the ordinary Emscripten buffer hook for our authenticated gzip
    # loader. The compiler WASM bytes remain exactly those of the pinned release.
    for before, after in [
        ('var wasmBinary;', 'var wasmBinary=Module.wasmBinary;'),
        ('async function getWasmBinary(binaryFile){', 'async function getWasmBinary(binaryFile){if(wasmBinary)return wasmBinary;'),
    ]:
        if glue.count(before) != 1:
            raise RuntimeError('Pinned Emscripten buffer adapter needs review')
        glue = glue.replace(before, after)
    (out/'lean.js').write_text(glue)
    (out/'SSC.lean').write_text(source)
    packs, parts, entries, size = [], [], [], 0
    def flush():
        nonlocal parts, entries, size
        if not parts:
            return
        name = f'library-{len(packs):03d}.pack'
        data = gzip.compress(b''.join(parts), compresslevel=9, mtime=0)
        (out/name).write_bytes(data)
        packs.append({'file': name, 'sha256': digest(out/name), 'compressedBytes': len(data),
                      'bytes': size, 'entries': entries})
        parts, entries, size = [], [], 0
    for file in sorted((kit/'lib/lean').rglob('*')):
        if not file.is_file() or not str(file).endswith(('.olean', '.ir', '.ir.sig')):
            continue
        rel = file.relative_to(kit/'lib/lean').as_posix()
        if rel.split('/')[0].split('.')[0] not in ['Init', 'Std', 'Lean', 'SSC']:
            continue
        data = file.read_bytes()
        if parts and size + len(data) > 16*1024*1024:
            flush()
        entries.append({'path': rel, 'offset': size, 'bytes': len(data)})
        parts.append(data); size += len(data)
    flush()
    manifest = {'schema': 1, 'leanCommit': LEAN_COMMIT, 'leanVersion': '4.33.0-pre',
                'sourceSha256': hashlib.sha256(source.encode()).hexdigest(), 'packs': packs,
                'runtime': {name: {'sha256': digest(out/name), 'bytes': (out/name).stat().st_size}
                            for name in ['lean.js', 'lean.wasm.pack', 'study.mjs', 'study.wasm', 'study-runtime.txt']}}
    manifest['version'] = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()[:24]
    manifest['downloadBytes'] = sum(p['compressedBytes'] for p in packs) + sum(p['bytes'] for p in manifest['runtime'].values())
    (out/'manifest.json').write_text(json.dumps(manifest, separators=(',', ':')))
    shutil.copyfile(ROOT/'packaging/lean-wasm/UPSTREAM-LICENSE', out/'UPSTREAM-LICENSE')
    print(f'Browser compiler: {len(packs)} packs, {manifest["downloadBytes"]/1048576:.1f} MiB compressed', flush=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path(tempfile.gettempdir())/'ssc-lean-wasm-cache')
    parser.add_argument('--output', type=Path, default=ROOT/'web/lean-wasm')
    parser.add_argument('--kit', type=Path, help='Use an already authenticated/extracted build kit')
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    kit = (args.kit or args.cache/'kit').resolve()
    if not args.kit:
        marker = kit/'verified-release.json'
        if not marker.exists() or json.loads(marker.read_text()) != {k: list(v) for k, v in RELEASES.items()}:
            unpack(download('runtime', args.cache), kit, True)
            unpack(download('static', args.cache), kit, False)
            marker.write_text(json.dumps(RELEASES))
    header = (kit/'lib/lean/Init.olean').read_bytes()[:160]
    if LEAN_COMMIT.encode() not in header:
        raise RuntimeError('Lean module/toolchain version mismatch')
    source = browser_source()
    build = args.cache/'build'
    build.mkdir(exist_ok=True)
    (build/'SSC.lean').write_text(source)
    print('Compiling the statistical module with Lean WebAssembly', flush=True)
    with (build/'compiler.log').open('w') as log:
        result = subprocess.run(node_command() + [str(ROOT/'packaging/lean-wasm/lean-node.cjs'), str(kit), str(build),
            '-R', '/work', '-o', '/work/SSC.olean', '/work/SSC.lean'], stdout=log, stderr=subprocess.STDOUT, timeout=300)
    if result.returncode != 0 or not (build/'SSC.olean').is_file():
        raise RuntimeError((build/'compiler.log').read_text()[-12000:])
    shutil.copyfile(build/'SSC.olean', kit/'lib/lean/SSC.olean')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    print('Compiling the C++ numerical backend to WebAssembly', flush=True)
    subprocess.run(['em++', '-std=c++17', '-O2', '-fwasm-exceptions', '-I'+str(ROOT/'core/src'),
        str(ROOT/'packaging/lean-wasm/study-runtime.cpp'), str(ROOT/'core/src/util/json.cpp'),
        str(ROOT/'core/src/util/numtext.cpp'), '-o', str(out/'study.mjs'),
        '-sMODULARIZE=1', '-sEXPORT_ES6=1', '-sENVIRONMENT=worker,node', '-sALLOW_MEMORY_GROWTH=1',
        '-sEXPORTED_FUNCTIONS=["_ssc_evaluate_study","_malloc","_free"]',
        '-sEXPORTED_RUNTIME_METHODS=["cwrap"]'], check=True)
    shutil.copyfile(ROOT/'packaging/lean-wasm/study-runtime.hpp', out/'study-runtime.txt')
    package(kit, out, source)

if __name__ == '__main__':
    main()
