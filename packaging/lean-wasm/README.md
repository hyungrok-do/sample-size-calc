# Static browser compiler

This is the default compiler used by `web/lean.html`. It requires only static
hosting; the Docker/native deployment is optional.

With the Emscripten 6.0.13 SDK on PATH, run `python scripts/build_lean_wasm.py`. The script
verifies immutable public runtime and library archives, builds the statistical
module using the actual WASM compiler, builds the shared C++ numerical runtime,
and packages all required files in `web/lean-wasm/`. `web/assemble.sh` includes
that directory, and the Pages workflow tests the assembled application before
publishing it.

The build uses the SDK's configured Node runtime. An older system Node can
reject Lean's export count before initialization, so it is not used implicitly.

`lean-node.cjs` is adapted from cauli/lean4-wasm-in-browser (Apache-2.0). The
pinned runtime and core library target wasm32. Ordinary native 64-bit `.olean`
files are not substituted for those artifacts. The statistical source adapter
uses the same model/elaborator definitions as the native project, adapts integer
quotation to the pinned core library, and separates the native Mathlib theorem
audit from the browser's per-study arithmetic certificates.

`study-runtime.hpp` is used both by the WASM numerical module and the standalone
C++ download. Statistical assumptions and floating-point power remain separate
from Lean's exact rational certificates. Full details and the browser test are
in [the compiler documentation](../../docs/lean-compiler.md).
