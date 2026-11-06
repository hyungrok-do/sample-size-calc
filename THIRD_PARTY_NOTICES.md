# Third-party software

Sample Size Calculator's own code is in this repository, under the MIT License
([`LICENSE`](LICENSE)). The C++ engine uses only the C++ standard library of
the compiler that builds it.

## The web page's engine

The engine files of the web page (`ssc_web.mjs`, `ssc_web.wasm` and their
multi-thread versions `ssc_web_mt.*`, built by `web/build.sh`) are compiled
with [Emscripten](https://emscripten.org) and contain parts of its runtime and
system libraries:

| Component | License |
|---|---|
| Emscripten's JavaScript runtime and support code | MIT or University of Illinois/NCSA Open Source License |
| musl libc | MIT |
| LLVM's libc++, libc++abi, libunwind and compiler-rt | Apache License 2.0 with LLVM Exceptions |
| dlmalloc (the memory allocator) | Public domain (CC0) |

The full license texts are in the Emscripten SDK
(`emscripten/LICENSE`, `emscripten/system/lib/*/LICENSE*`) and at
<https://github.com/emscripten-core/emscripten>.

## Lean browser compiler

The browser compiler build downloads the Lean 4 WebAssembly runtime from
[cauli/lean4-wasm-in-browser](https://github.com/cauli/lean4-wasm-in-browser),
built from [cauli/lean4](https://github.com/cauli/lean4) at
`62b6a2291302d4bbeace37642a066b7510d0145c`, a fork of the
[Lean 4 theorem prover](https://github.com/leanprover/lean4). Lean and the browser
project use Apache-2.0. Exact release URLs and SHA-256 hashes are recorded in
`scripts/build_lean_wasm.py`. The adapted Node build driver and runtime
initialization acknowledge that project; its license is included in
`packaging/lean-wasm/UPSTREAM-LICENSE` and the deployed compiler directory.

The compact asset loader uses [brotli-dec-wasm](https://github.com/ustclug-dev/brotli-dec-wasm)
2.3.2 (MIT OR Apache-2.0), a decoder-only WebAssembly wrapper around the Rust
`brotli-decompressor` implementation. It is distributed here under Apache-2.0.
Its release archive is SHA-256 pinned in the build script. The deployed
`lean-wasm/brotli/LICENSE-Apache.txt` contains the full license; the original MIT
license is also included as `LICENSE-MIT.txt`. This temporary decoder
worker performs lossless asset decompression, and is terminated before Lean
imports its modules. It does not replace any part of the Lean compiler or kernel.
