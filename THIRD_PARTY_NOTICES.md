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
