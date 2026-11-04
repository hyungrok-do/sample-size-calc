# Static browser compiler

`web/lean.html` runs ordinary Lean in browser WebAssembly. Build with Node 24+
using `python scripts/build_lean_wasm.py`. If Emscripten is on PATH, the builder
uses its configured Node runtime. The Lean compiler tab itself requires no
Emscripten or separate C++ numerical module.

The builder verifies pinned public runtime and library archives, compiles
`lean/browser/SSC.lean` with the actual WASM compiler, and packages matching
wasm32 artifacts in `web/lean-wasm/`. Manifest schema 2 contains the Lean runtime,
core libraries and proof reporter. The old statistical module is not imported.

`lean-node.cjs` is adapted from cauli/lean4-wasm-in-browser (Apache-2.0). Native
64-bit `.olean` files must not be substituted. `SSC.lean` contains a generic,
proved finite search and observes ordinary Lean declarations/tactics; it has
no statistical model dispatch or formula templates.

`web/assemble.sh` includes the package. The Pages workflow tests the actual
assembled compiler before publishing. The calculator and DAG tabs still build
their independent engines with Emscripten 6.0.13.

See [compiler scope and verification](../../docs/lean-compiler.md).
