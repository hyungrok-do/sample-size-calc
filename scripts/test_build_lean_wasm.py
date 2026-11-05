"""Cache invalidation and failure atomicity, without invoking a large compiler."""
import importlib.util
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location('builder', Path(__file__).with_name('build_lean_wasm.py'))
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


class BuildCache(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.kit, self.build, self.cache = [self.root/p for p in ['kit', 'build', 'cache']]
        self.build.mkdir()
        self.calls = 0

    def compiler(self, command, **kwargs):
        self.calls += 1
        output = self.build/command[command.index('-o')+1].removeprefix('/work/')
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(builder.LEAN_COMMIT.encode()+str(self.calls).encode())
        return SimpleNamespace(returncode=0)

    def run_module(self, source='def value := 1', inputs=None):
        return builder.compile_module(self.kit, self.build, self.cache, 'SSC.lean', source, inputs or ['toolchain'])

    @patch.object(builder, 'node_command', return_value=['node'])
    def test_reuse_and_all_input_changes(self, _):
        with patch.object(builder.subprocess, 'run', side_effect=self.compiler):
            first = self.run_module()
            self.assertEqual(self.run_module(), first)
            self.assertEqual(self.calls, 1)
            self.run_module(source='def value := 2')
            self.run_module(inputs=['different-toolchain'])
            self.run_module(inputs=['toolchain', 'changed-dependency-output'])
            self.assertEqual(self.calls, 4)

    @patch.object(builder, 'node_command', return_value=['node'])
    def test_corruption_rebuilds(self, _):
        with patch.object(builder.subprocess, 'run', side_effect=self.compiler):
            self.run_module()
            next(p for p in self.cache.iterdir() if not p.suffix).write_bytes(b'corrupt')
            self.run_module()
            self.assertEqual(self.calls, 2)

    @patch.object(builder, 'node_command', return_value=['node'])
    def test_failure_cannot_reuse_stale_output(self, _):
        (self.build/'SSC.olean').write_bytes(builder.LEAN_COMMIT.encode())
        with patch.object(builder.subprocess, 'run', return_value=SimpleNamespace(returncode=0)):
            with self.assertRaises(RuntimeError):
                self.run_module()
        self.assertFalse((self.kit/'lib/lean/SSC.olean').exists())

    @patch.object(builder, 'node_command', return_value=['node'])
    def test_error_log_overrides_success_exit(self, _):
        def failed(command, **kwargs):
            result = self.compiler(command, **kwargs)
            kwargs['stdout'].write('Reporting complete, hasErrors = true\n')
            return result
        with patch.object(builder.subprocess, 'run', side_effect=failed):
            with self.assertRaises(RuntimeError):
                self.run_module()
        self.assertFalse((self.build/'SSC.olean').exists())
        self.assertFalse(self.cache.exists())


class StablePacks(unittest.TestCase):
    def test_ssc_growth_preserves_standard_library_and_pack_contents(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            kit, out = root/'kit', root/'out'
            (kit/'bin').mkdir(parents=True)
            (kit/'bin/lean.wasm').write_bytes(b'compiler fixture')
            (kit/'bin/lean.js').write_text('var wasmBinary;async function getWasmBinary(binaryFile){}')
            library = kit/'lib/lean'
            library.mkdir(parents=True)
            for group in ['Init', 'Lean', 'SSC', 'Std']:
                (library/f'{group}.olean').write_bytes((group * 100).encode())
            sources = {'SSC.lean': 'import Lean'}
            builder.package(kit, out, sources)
            original = json.loads((out/'manifest.json').read_text())
            # Cross a pack-size boundary in SSC and move between the old/new
            # naming schemes. Neither event may invalidate unrelated packs.
            (library/'SSC').mkdir()
            (library/'SSC/Large.olean').write_bytes(b'x' * (16*1024*1024 + 1))
            (out/'library-999.pack').write_bytes(b'obsolete')
            builder.package(kit, out, sources)
            updated = json.loads((out/'manifest.json').read_text())
            before = {p['file']: p['sha256'] for p in original['packs'] if not p['file'].startswith('library-ssc-')}
            after = {p['file']: p['sha256'] for p in updated['packs'] if not p['file'].startswith('library-ssc-')}
            self.assertEqual(before, after)
            self.assertFalse((out/'library-999.pack').exists())
            recovered = {}
            for pack in updated['packs']:
                data = gzip.decompress((out/pack['file']).read_bytes())
                self.assertEqual(len(data), pack['bytes'])
                for entry in pack['entries']:
                    start = entry['offset']
                    recovered[entry['path']] = data[start:start+entry['bytes']]
            self.assertEqual(recovered, {p.relative_to(library).as_posix(): p.read_bytes()
                                         for p in library.rglob('*.olean')})


if __name__ == '__main__':
    unittest.main()
