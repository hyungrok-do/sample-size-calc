#!/usr/bin/env sh
# Assembles the site as published: web/assemble.sh <new site dir> <version>.  The Pages workflow publishes it and
# CI tests it with web/tests/site.mjs.  It holds the page, every script at the top of web/, the two engines (built
# by web/build.sh), the license and the notices, with the version stamped into the addresses (web/stamp.sh).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
repo=$(dirname "$here")
site=$1
version=$2
[ -f "$here/engine/ssc_web.wasm" ] && [ -f "$here/engine/ssc_web_mt.wasm" ] ||
  { echo "web/engine has no engines: run web/build.sh first" >&2; exit 1; }
if [ -e "$site" ] && [ -n "$(ls -A "$site")" ]; then echo "$site exists and is not empty" >&2; exit 1; fi
mkdir -p "$site/engine"
cp "$here/index.html" "$here"/*.js "$site/"
cp "$here/engine/ssc_web.mjs" "$here/engine/ssc_web.wasm" "$here/engine/ssc_web_mt.mjs" "$here/engine/ssc_web_mt.wasm" "$site/engine/"
cp "$repo/LICENSE" "$repo/THIRD_PARTY_NOTICES.md" "$site/"
"$here/stamp.sh" "$site" "$version"
