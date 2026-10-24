// Whether the engines in web/engine were built from the sources beside them.  web/build.sh writes the commit it
// built from and a hash of the changes from it to the engine's files, and of the files there that git does not track,
// to engine/source.txt; the tests compare them with the sources here and warn when they differ, since they then check
// an engine built from other code.
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const web = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const repo = resolve(web, '..');
const SOURCES = ['CMakeLists.txt', 'core/src', 'core/include'];  // as in web/build.sh

// Why the engines may not match the sources ('' when they do).
export function engineMismatch() {
  let commit, changes;
  try {
    [commit, changes] = readFileSync(resolve(web, 'engine/source.txt'), 'utf8').trim().split(' ');
  } catch {
    return 'web/engine does not say which sources its engines were built from (web/build.sh writes engine/source.txt)';
  }
  let now;
  try {
    const git = (args, input) => execFileSync('git', args,
      { cwd: repo, input, stdio: [input ? 'pipe' : 'ignore', 'pipe', 'ignore'], maxBuffer: 1 << 30 });
    const diff = git(['diff', '--no-ext-diff', '--no-color', '--binary', commit, '--', ...SOURCES]);
    // CMake builds every core/src/*.cpp, tracked or not
    const untracked = git(['ls-files', '-z', '--others', '--exclude-standard', '--', ...SOURCES]);
    const names = untracked.toString().split('\0').filter(Boolean);
    const contents = names.length ? git(['hash-object', '--', ...names]) : Buffer.alloc(0);
    now = git(['hash-object', '--stdin'], Buffer.concat([diff, untracked, contents])).toString().trim();
  } catch {
    return `the engines in web/engine were built from commit ${commit.slice(0, 12)}, which git cannot compare with these sources`;
  }
  return now === changes ? ''
    : `the engines in web/engine were built from other sources than these (from commit ${commit.slice(0, 12)}): `
      + 'rebuild them with web/build.sh';
}

// A warning on the test's output, not a failure: the engines take a while to build, and need the Emscripten SDK.
export function warnIfOtherSources() {
  const why = engineMismatch();
  if (why) console.warn(`warning: ${why}`);
}
