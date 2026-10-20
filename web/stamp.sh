#!/usr/bin/env sh
# Puts a version into the addresses the published page loads its files from: web/stamp.sh <site dir> <version>.
# The source uses ?v=dev (and the cache name ssc-dev); a browser keeps files for a while, and a new address for
# each version keeps it from running a new page with an older script.
set -eu
site=$1
version=$2
for file in index.html app.js engine.js worker.js sw.js; do
  grep -q 'v=dev' "$site/$file" || { echo "$file: no ?v=dev to stamp" >&2; exit 1; }
  sed -i "s/v=dev/v=$version/g; s/'ssc-dev'/'ssc-$version'/" "$site/$file"
done
if grep -n 'v=dev\|ssc-dev' "$site"/*.html "$site"/*.js; then
  echo "left unstamped (above)" >&2
  exit 1
fi
echo "stamped version $version"
