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
import re
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import time

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

def browser_sources():
    # Ordinary modules in dependency order, compiled by the shipped compiler.
    # No model adapters or source-to-source substitutions.
    return {name: (ROOT/'lean/browser'/name).read_text() for name in
            ['SSC/Numeric.lean', 'SSC/Probability.lean', 'SSC/Calculation.lean',
             'SSC/Finite.lean', 'SSC/Derivation.lean', 'SSC.lean']}


def node_command():
    # emsdk_env.sh adds emcc to PATH, but can leave an older system Node first.
    # Use the SDK's configured engine, as emcc does. Older V8 versions reject
    # the genuine Lean runtime's 104769 exports before running any Lean code.
    if not shutil.which('em-config'):
        return ['node']  # Node 24+ is required by the pinned runtime's export count.
    command = ast.literal_eval(subprocess.check_output(['em-config', 'NODE_JS'], text=True).strip())
    if isinstance(command, str):
        command = [command]
    if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
        raise RuntimeError('Emscripten NODE_JS must name a JavaScript runtime')
    return command

def toolchain_digest(kit):
    """Fingerprint actual compiler inputs, including a caller-supplied kit."""
    files = [kit/'bin/lean.js', kit/'bin/lean.wasm', ROOT/'packaging/lean-wasm/lean-node.cjs',
             Path(__file__).resolve()]
    files += sorted(p for p in (kit/'lib/lean').rglob('*') if p.is_file()
                    and p.suffix in ['.olean', '.ir', '.sig']
                    and p.relative_to(kit/'lib/lean').parts[0].split('.')[0] in ['Init', 'Std', 'Lean'])
    return hashlib.sha256(json.dumps([(str(p.relative_to(kit)) if p.is_relative_to(kit) else p.name,
                                       digest(p)) for p in files]).encode()).hexdigest()

def compile_module(kit, build, cache, name, source, inputs):
    # Include prior outputs, not just source hashes. This conservatively
    # invalidates downstream modules even after a toolchain or dependency change.
    key = hashlib.sha256(json.dumps([1, inputs, name, source]).encode()).hexdigest()
    cached = cache/key
    artifact = Path(name).with_suffix('.olean')
    destination = kit/'lib/lean'/artifact
    destination.parent.mkdir(parents=True, exist_ok=True)
    metadata = cached.with_suffix('.json')
    if cached.is_file() and metadata.is_file():
        try:
            expected = json.loads(metadata.read_text())['sha256']
            if digest(cached) == expected:
                shutil.copyfile(cached, destination)
                print(f'Reusing checked build artifact {name}', flush=True)
                return expected
        except (ValueError, KeyError):
            pass
    path = build/name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(source)
    output = build/artifact
    # A failed compiler invocation must never reuse yesterday's output.
    output.unlink(missing_ok=True)
    started = time.monotonic()
    print(f'Compiling {name} with Lean WebAssembly', flush=True)
    with (build/'compiler.log').open('w') as log:
        result = subprocess.run(node_command() + [str(ROOT/'packaging/lean-wasm/lean-node.cjs'), str(kit), str(build),
            '-R', '/work', '-o', '/work/'+artifact.as_posix(), '/work/'+name],
            stdout=log, stderr=subprocess.STDOUT, timeout=300)
    log = (build/'compiler.log').read_text()
    if result.returncode != 0 or not output.is_file() or re.search(r'hasErrors\s*=\s*true', log):
        output.unlink(missing_ok=True)
        raise RuntimeError(log[-12000:])
    if LEAN_COMMIT.encode() not in output.read_bytes()[:160]:
        raise RuntimeError(f'Compiler produced an incompatible artifact for {name}')
    sha = digest(output)
    cache.mkdir(parents=True, exist_ok=True)
    temporary = cached.with_suffix('.tmp')
    shutil.copyfile(output, temporary)
    temporary.replace(cached)
    metadata.write_text(json.dumps({'sha256': sha}))
    shutil.copyfile(output, destination)
    print(f'Compiled {name} in {time.monotonic()-started:.1f}s', flush=True)
    return sha

def package(kit, out, sources):
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
    for name, source in sources.items():
        (out/name).parent.mkdir(parents=True, exist_ok=True)
        (out/name).write_text(source)
    packs, parts, entries, size = [], [], [], 0
    group, group_index = None, 0
    def flush():
        nonlocal parts, entries, size, group_index
        if not parts:
            return
        # Changing SSC must not move every subsequent Std byte into a different
        # pack. Stable namespace boundaries and names preserve content-cache hits.
        name = f'library-{group.lower()}-{group_index:03d}.pack'
        group_index += 1
        data = gzip.compress(b''.join(parts), compresslevel=9, mtime=0)
        (out/name).write_bytes(data)
        packs.append({'file': name, 'sha256': digest(out/name), 'compressedBytes': len(data),
                      'bytes': size, 'entries': entries})
        parts, entries, size = [], [], 0
    for file in sorted((kit/'lib/lean').rglob('*')):
        if not file.is_file() or not str(file).endswith(('.olean', '.ir', '.ir.sig')):
            continue
        rel = file.relative_to(kit/'lib/lean').as_posix()
        next_group = rel.split('/')[0].split('.')[0]
        if next_group not in ['Init', 'Std', 'Lean', 'SSC']:
            continue
        if next_group != group:
            flush()
            group, group_index = next_group, 0
        data = file.read_bytes()
        if parts and size + len(data) > 16*1024*1024:
            flush()
        entries.append({'path': rel, 'offset': size, 'bytes': len(data)})
        parts.append(data); size += len(data)
    flush()
    manifest = {'schema': 2, 'leanCommit': LEAN_COMMIT, 'leanVersion': '4.33.0-pre',
                'sourceSha256': hashlib.sha256(sources['SSC.lean'].encode()).hexdigest(),
                'librarySources': {name: hashlib.sha256(source.encode()).hexdigest()
                                   for name, source in sources.items()}, 'packs': packs,
                'runtime': {name: {'sha256': digest(out/name), 'bytes': (out/name).stat().st_size}
                            for name in ['lean.js', 'lean.wasm.pack']}}
    manifest['version'] = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()[:24]
    manifest['downloadBytes'] = sum(p['compressedBytes'] for p in packs) + sum(p['bytes'] for p in manifest['runtime'].values())
    (out/'manifest.json').write_text(json.dumps(manifest, separators=(',', ':')))
    active = {pack['file'] for pack in packs}
    for obsolete in out.glob('library-*.pack'):
        if obsolete.name not in active:
            obsolete.unlink()
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
    sources = browser_sources()
    build = args.cache/'build'
    build.mkdir(exist_ok=True)
    inputs = [toolchain_digest(kit)]
    for name, source in sources.items():
        inputs.append(compile_module(kit, build, args.cache/'modules', name, source, inputs))
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    for obsolete in ['study.mjs', 'study.wasm', 'study-runtime.txt']:
        (out/obsolete).unlink(missing_ok=True)
    package(kit, out, sources)

if __name__ == '__main__':
    main()
