#!/usr/bin/env sh
# Builds the engine for the web page into web/engine/:
#   ssc_web.mjs + ssc_web.wasm        one thread; works from any web server
#   ssc_web_mt.mjs + ssc_web_mt.wasm  every core; needs a cross-origin isolated page (see serve.py)
# Needs the Emscripten SDK on PATH: source <emsdk>/emsdk_env.sh first.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
repo=$(dirname "$here")
command -v emcmake >/dev/null || { echo "emcmake not found: source <emsdk>/emsdk_env.sh first" >&2; exit 1; }
mkdir -p "$here/engine"
for threads in OFF ON; do
  dir="$repo/build-wasm-$(echo "$threads" | tr 'A-Z' 'a-z')"
  emcmake cmake -S "$repo" -B "$dir" -DCMAKE_BUILD_TYPE=Release -DSSC_BUILD_TESTS=OFF -DSSC_WASM_THREADS="$threads"
  cmake --build "$dir" --target ssc_web --parallel
  cp "$dir"/ssc_web*.mjs "$dir"/ssc_web*.wasm "$here/engine/"
done
ls -l "$here/engine"
