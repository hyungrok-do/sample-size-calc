#!/usr/bin/env sh
# Puts a version into the addresses the published page loads its files from: web/stamp.sh <site dir> <version>.
# The source uses ?v=dev (and the cache name ssc-dev); a browser keeps files for a while, and a new address for
# each version keeps it from running a new page with an older script.  The page's <style> holds no ?v=dev and is left
# as it is: the page's Content-Security-Policy names it by its hash (tests/site.mjs checks it after stamping).
set -eu
site=$1
version=$2
for file in index.html app.js engine.js worker.js sw.js guide.js; do
  grep -q 'v=dev' "$site/$file" || { echo "$file: no ?v=dev to stamp" >&2; exit 1; }
  sed -i "s/v=dev/v=$version/g; s/'ssc-dev'/'ssc-$version'/" "$site/$file"
done
# exactly v=dev, so a version such as dev-test is not taken for it
if grep -nE "v=dev([^[:alnum:]_.-]|$)|'ssc-dev'" "$site"/*.html "$site"/*.js; then
  echo "left unstamped (above)" >&2
  exit 1
fi
echo "stamped version $version"
