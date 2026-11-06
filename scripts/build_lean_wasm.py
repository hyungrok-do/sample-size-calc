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
import sys
import tarfile
import tempfile
import time
import zlib

ROOT = Path(__file__).resolve().parents[1]
LEAN_COMMIT = '62b6a2291302d4bbeace37642a066b7510d0145c'
RELEASES = {
    'runtime': ('https://github.com/cauli/lean4-wasm-in-browser/releases/download/runtime-62b6a22-editor1/lean-runtime-fixture.tar.gz',
                'a66af27e859464db9849df10ea3bd7aa3b5ee824fbcc7c330fb8cf5931cc5527'),
    'static': ('https://github.com/cauli/lean4-wasm-in-browser/releases/download/runtime-62b6a22-compact1/pages-assets-62b6a22-compact1.tar.gz',
               'e0c1ae4d40fd5c4da18d4a38e254e8adef603b855757daff103659a4b80c9a39'),
    'brotli': ('https://registry.npmjs.org/brotli-dec-wasm/-/brotli-dec-wasm-2.3.2.tgz',
               '5006ea904407e2f44ad6616100de196c39af591b772f764ba66fdd720b5d2ab5'),
}
BROTLI_FILES = {'module': 'brotli_dec_wasm.js', 'wasm': 'brotli_dec_wasm_bg.wasm'}
BROTLI_LICENSES = ['LICENSE-Apache.txt', 'LICENSE-MIT.txt']
OLD_BROTLI_FILES = ['brotli_wasm.js', 'brotli_wasm_bg.wasm', 'LICENSE']

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


def unpack_brotli(archive, kit):
    """Extract only the pinned browser decoder and its license."""
    expected = {f'package/pkg/{name}': name for name in BROTLI_FILES.values()}
    expected.update({f'package/{name}': name for name in BROTLI_LICENSES})
    found = set()
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            name = member.name.removeprefix('./')
            if name not in expected:
                continue
            if not member.isfile() or name in found:
                raise RuntimeError(f'Unexpected decoder archive entry: {name}')
            destination = kit/'brotli'/expected[name]
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as source, destination.open('wb') as output:
                shutil.copyfileobj(source, output)
            found.add(name)
    if found != set(expected):
        raise RuntimeError('Pinned Brotli decoder archive is incomplete')
    for name in OLD_BROTLI_FILES:
        (kit/'brotli'/name).unlink(missing_ok=True)

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

def _atomic_write(path, data):
    """Publish complete cache records even if a build is interrupted."""
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=path.name+'.', delete=False) as f:
        temporary = Path(f.name)
        try:
            f.write(data)
            f.flush()
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)


def compressed_asset(data, destination, cache, codec, compressor_identity):
    """Lossless compression cached by exact bytes, options and compressor version.

    The cache lives outside the published directory. Every compressed cache hit
    is hash-verified; a corrupt/missing record is rebuilt. Existing gzip assets
    can seed a fresh cache only after their complete decoded bytes are verified.
    """
    raw_sha = hashlib.sha256(data).hexdigest()
    inputs = [1, codec, compressor_identity, raw_sha]
    key = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
    cache.mkdir(parents=True, exist_ok=True)
    cached, metadata = cache/key, cache/(key+'.json')
    compressed = None
    if cached.is_file() and metadata.is_file():
        try:
            record = json.loads(metadata.read_text())
            candidate = cached.read_bytes()
            if record['inputs'] == inputs and hashlib.sha256(candidate).hexdigest() == record['sha256']:
                compressed = candidate
        except (ValueError, KeyError, OSError, TypeError):
            pass
    if compressed is None:
        # Preserve previous deterministic gzip artifacts without spending minutes
        # recompressing unchanged libraries on the first cache-enabled build.
        if codec == 'gzip' and destination.is_file():
            try:
                candidate = destination.read_bytes()
                if hashlib.sha256(gzip.decompress(candidate)).hexdigest() == raw_sha:
                    compressed = candidate
            except (OSError, EOFError, zlib.error):
                pass
        if compressed is None:
            if codec == 'gzip':
                compressed = gzip.compress(data, compresslevel=9, mtime=0)
            elif codec == 'brotli':
                script = '''const fs=require('node:fs'),z=require('node:zlib');
const data=fs.readFileSync(0);
process.stdout.write(z.brotliCompressSync(data,{params:{
  [z.constants.BROTLI_PARAM_QUALITY]:6,
  [z.constants.BROTLI_PARAM_SIZE_HINT]:data.length
}}));'''
                compressed = subprocess.check_output(node_command()+['-e', script], input=data)
            else:
                raise ValueError(f'Unsupported compression codec: {codec}')
        _atomic_write(cached, compressed)
        _atomic_write(metadata, json.dumps({'inputs': inputs,
            'sha256': hashlib.sha256(compressed).hexdigest()}).encode())
    destination.write_bytes(compressed)
    return {'file': destination.name, 'sha256': hashlib.sha256(compressed).hexdigest(),
            'bytes': len(compressed)}


def package(kit, out, sources, compression_cache=None):
    out.mkdir(parents=True, exist_ok=True)
    compression_cache = compression_cache or kit.parent/'compressed'
    # The existing build already requires Node 24+. Its built-in Brotli codec
    # avoids a Python package dependency or an external compression executable.
    brotli_identity = json.loads(subprocess.check_output(node_command()+['-p',
        'JSON.stringify([process.versions.node,process.versions.brotli])'], text=True))
    identities = {'gzip': [9, 0, list(sys.version_info[:3]), zlib.ZLIB_RUNTIME_VERSION],
                  'brotli': [6, brotli_identity]}

    def compress(data, filename):
        gzip_info = compressed_asset(data, out/filename, compression_cache, 'gzip', identities['gzip'])
        brotli_name = filename.removesuffix('.pack')+'.br.pack'
        gzip_info['brotli'] = compressed_asset(data, out/brotli_name, compression_cache,
                                               'brotli', identities['brotli'])
        return gzip_info

    # Compressed opaque extensions avoid a static host decompressing them twice.
    wasm = compress((kit/'bin/lean.wasm').read_bytes(), 'lean.wasm.pack')
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
        data = b''.join(parts)
        compressed = compress(data, name)
        packs.append({'file': name, 'sha256': compressed['sha256'],
                      'compressedBytes': compressed['bytes'], 'brotli': compressed['brotli'],
                      'uncompressedSha256': hashlib.sha256(data).hexdigest(),
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
                'runtime': {'lean.js': {'sha256': digest(out/'lean.js'),
                                       'bytes': (out/'lean.js').stat().st_size},
                            'lean.wasm.pack': {k: v for k, v in wasm.items() if k != 'file'}}}
    manifest['runtime']['lean.wasm.pack'].update({
        'uncompressedSha256': digest(kit/'bin/lean.wasm'),
        'uncompressedBytes': (kit/'bin/lean.wasm').stat().st_size})
    if (kit/'brotli').is_dir():
        decoder = {}
        for key, name in BROTLI_FILES.items():
            source = kit/'brotli'/name
            sha = digest(source)
            # The helper imports this URL after the parent has verified it.
            # An immutable module URL prevents a deployment from replacing its
            # contents between those two requests without relaxing script CSP.
            published_name = f'{Path(name).stem}-{sha}{Path(name).suffix}' if key == 'module' else name
            path = out/'brotli'/published_name
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, path)
            decoder[key] = {'file': path.relative_to(out).as_posix(),
                            'sha256': sha, 'bytes': path.stat().st_size}
        for name in BROTLI_LICENSES:
            shutil.copyfile(kit/'brotli'/name, out/'brotli'/name)
        active_decoder = set(BROTLI_LICENSES) | {Path(part['file']).name for part in decoder.values()}
        for obsolete in (out/'brotli').iterdir():
            if obsolete.is_file() and obsolete.name not in active_decoder:
                obsolete.unlink()
        manifest['decoder'] = decoder
    manifest['version'] = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()[:24]
    manifest['downloadBytes'] = sum(p['compressedBytes'] for p in packs) + sum(p['bytes'] for p in manifest['runtime'].values())
    manifest['brotliDownloadBytes'] = sum(p['brotli']['bytes'] for p in packs) + \
        wasm['brotli']['bytes'] + manifest['runtime']['lean.js']['bytes'] + \
        sum(p['bytes'] for p in manifest.get('decoder', {}).values())
    (out/'manifest.json').write_text(json.dumps(manifest, separators=(',', ':')))
    active = {name for pack in packs for name in [pack['file'], pack['brotli']['file']]}
    for obsolete in out.glob('library-*.pack'):
        if obsolete.name not in active:
            obsolete.unlink()
    shutil.copyfile(ROOT/'packaging/lean-wasm/UPSTREAM-LICENSE', out/'UPSTREAM-LICENSE')
    print(f'Browser compiler: {len(packs)} packs, '
          f'{manifest["downloadBytes"]/1048576:.1f} MiB gzip / '
          f'{manifest["brotliDownloadBytes"]/1048576:.1f} MiB Brotli', flush=True)

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
            unpack_brotli(download('brotli', args.cache), kit)
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
    package(kit, out, sources, args.cache/'compressed')

if __name__ == '__main__':
    main()
