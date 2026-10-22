# Web page

The calculator's C++ engine compiled to WebAssembly, with a page that runs it
in the browser. You can search the methods, fill in a form generated from each
method's description and calculate. Nothing you type leaves the browser, and
the server only hands out files.

The page opens on **Browse all methods**, the searchable catalog. **Guide me**,
an experimental feature, asks study-design questions using chips. Each
question selects one answer; the optional additional requirements select all
that apply. From the first answer, **Recommended analysis** lists every
analysis that still fits, narrowed by each further answer and by the
requirements chosen; any of them opens in the calculator with the settings its
remaining routes agree on. "I'm unsure" lists what is possible from that
question on. A completed route shows its recommendation at once, and
**Use this analysis** opens the generated calculator below with the selected
analysis settings. Numerical defaults are example values to review, not
assumptions inferred from the guide.

The guide covers all 71 procedures. Unsupported and incompatible selections do
not open a recommended calculator. Changing an earlier answer clears dependent
choices, retains additional requirements, and hides the previous guided form.
Requirements are checked together against the selected procedure; they do not
automatically reroute to another design.
The guide is curated planning guidance, not a general statistical model selector.

`guide-data.js` contains the questions, procedure routes, compatibility lists,
and explanations; `guide-state.js` handles the selections; `guide.js` renders
them. The engine remains authoritative for parameter names, allowed options,
visibility rules, and calculations. `selection.js` validates a recommendation
before the form applies it. Reapplying the same procedure keeps numerical
study inputs when their context (design, effect scale, distribution, etc.) is
unchanged, and resets analysis options to the schema defaults plus the explicit
recommendation. A changed context starts with example values again. Changes to the form invalidate previous results,
including results still being calculated.

Existing links such as `#means.two.t` still work. Selected options can travel
in the fragment, for example `#prop.two?test=fisher&method=exact`; they are
validated against the engine schema on arrival. These links contain analysis
options, not study inputs or the guide's full answer history.

It is published at https://hyungrok-do.github.io/sample-size-calc/ by
`.github/workflows/pages.yml`. To run it on your own computer:

```bash
source /path/to/emsdk/emsdk_env.sh   # the Emscripten SDK: https://emscripten.org
web/build.sh                          # builds web/engine/ssc_web*.mjs and .wasm
python3 web/serve.py                  # then open http://localhost:8000
```

## Two engines

| Engine | Threads | Needs |
|---|---|---|
| `ssc_web` | 1 | any web server, including static hosting |
| `ssc_web_mt` | every core | a cross-origin isolated page: `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp` |

The page uses the threaded engine when it is cross-origin isolated, and the
other one otherwise. It also falls back to the one-thread engine if the
threaded one fails or has not started 20 s after compiling. Normally that takes
a tenth of a second. `?engine=single` forces the one-thread engine.

`serve.py` sends the two headers. A static host such as GitHub Pages cannot,
so the page's service worker (`sw.js`) adds them to every response it passes
on. The first visit therefore uses one thread, and every later visit uses every
core.

## Offline

The service worker keeps a copy of each file the page loads, the one-thread
engine included. The network always comes first, so a new version is used as
soon as it is published. Without a network the page starts from the copies,
and Stop restarts the engine from them. The multi-core engine is kept from the
first visit that uses it.

## Versions

Browsers keep a static host's files for a while (GitHub Pages allows 10
minutes), so a reload can pair a new page with an older script. To prevent
that, the page asks for every other file with the version in its address. The
source uses `?v=dev`, and `web/stamp.sh` replaces it with the commit when the
Pages workflow publishes the site. The multi-core engine starts its helper
threads from the same versioned address (`web/build.sh`), and the service
worker always checks the page itself with the server.

## Responsiveness and memory

- The engine runs in a Web Worker, so the page stays responsive during a long
  simulation. Searching and opening methods keep working while it runs: the
  page loads every method's description at the start, and it searches names
  and keywords itself while the engine is busy. Stop ends the worker and
  starts a fresh one, which takes about a tenth of a second.
- An answer can have up to 20,000 scenarios. The page shows the first 1,000,
  drawing a hundred rows per frame. With 20,000 scenarios, the page's longest
  pause is about 0.1 s, against 8.5 s when every row was drawn at once.
- The one-thread engine uses about 21 MB, and 37 MB for the largest simulated
  study it allows. Simulations keep only a block of studies at a time. The
  multi-core engine starts with 128 MB of shared memory, more than the largest
  answer needs (20,000 scenarios: 100 MB), because Safari's engine (WebKit)
  stopped it the first time its shared memory had to grow. WebAssembly memory
  never shrinks, so after an answer that grows the engine beyond its start
  (and beyond 64 MB), the page replaces the engine once it is idle and online.
- If the engine stops during a calculation (it traps, for example when out of
  memory), the page starts a fresh one at once. After the multi-core engine,
  it uses the one-thread engine for the rest of the visit and runs the
  calculation again; otherwise it explains what happened.
- No leaks were found:
  - The engine's allocated memory is identical after every round of all
    methods, the error cases and requests the engine refuses.
  - The C++ tests and the C interface run clean under AddressSanitizer and
    LeakSanitizer, including simulations cancelled while running.
  - The page holds steady over hundreds of calculations, method switches and
    searches. That needed one workaround: Chromium's style engine can keep a
    reference to an option of a removed `<select>`, which kept every old
    form alive. Emptying the selects before replacing a form fixes it.
- After a Stop, the threaded engine's helper threads take about 10 s to be
  freed by the browser.

## Measurements

These were measured on a 4-core Linux machine, comparing the page in headless
Chromium with the native library.

| | Browser, 1 thread | Browser, 4 threads | Native, 1 thread | Native, 4 threads |
|---|---|---|---|---|
| Every method and solve-for choice with its defaults (152 calculations, power curves included): median | 1.4 ms | 1.6 ms | 0.55 ms | 0.57 ms |
| … slowest | 101 ms | 46 ms | 79 ms | 31 ms |
| Win ratio, power, 1,000 simulated studies | 0.81 s | 0.29 s | 0.63 s | 0.17 s |
| Win ratio, sample size search, 1,000 studies per step | 11.1 s | 3.2 s | 8.7 s | 2.2 s |

The engine file (`ssc_web.wasm`) is 3.6 MB, or 1.0 MB with gzip, which
static hosts apply on their own. The page is ready about 0.1 s after the
download, once the browser has compiled the engine.

## Numbers

The page formats numbers as the command line does (`python/ssc/formatting.py`):
the same rounding, with exact ties going to the even digit, and the same
exponents. Over 56,000 values and every column format, the only difference is
that the page shows a negative zero as `0`.

## Tests

After building the web engines, run the frontend checks with Node 24:

```bash
npm ci --prefix web
cd web
npx playwright install chromium
npm test                      # the guide's routes and the number formatting (needs python3)
npm run test:browser          # the guide and calculator in a browser
npm run test:page             # the testers' fixes, responsiveness, the engine fallback
./assemble.sh /tmp/site dev-test && npm run test:site -- /tmp/site     # the site as published
```

`npm test` checks every guide route against the real WebAssembly catalog and
schemas: all-procedure coverage, the analyses listed after partial answers,
combined requirements, ancestor resets, unknown/conditional presets, and link
round trips. It also formats about 10,000 values in each column format on the
page and with the command line (`python/ssc/formatting.py`), and they must
agree. `test:browser` loads the actual page and engine, applies every
recommendation path, and checks a calculation, same-procedure option changes,
stale-result protection, keyboard use, direct links and history, mobile and dark
layout, and offline reload. `test:page` serves the page cross-origin isolated
(the multi-core engine) and checks each fix from the testers' reports, that a
20,000-scenario answer never stalls the page for long, and that the one-thread
engine takes over when the multi-core one is missing or its threads never start.

`web/assemble.sh` builds the site exactly as the Pages workflow publishes it.
`test:site` checks that every file the page can load was published with the
version in its address, then serves the site the way GitHub Pages does (a
subpath, no isolation headers, files kept for 10 minutes). It checks a first
visit, the guide, every core from the second visit, a start without a network,
and an update that must not mix the new page with older files. The Pages
workflow runs it before publishing.

CI runs all of these in Chromium, and `test:browser`, `test:page` and
`test:site` in Firefox and WebKit (Safari's engine) too. `SSC_BROWSER` selects
`chromium`, `firefox` or `webkit`. `PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH` can
select an installed Chromium binary, and `SSC_SCREENSHOTS` optionally saves
layout screenshots to a directory.

The engine's own tests run under Node with the threaded build (some of them
test the thread pool). They need a build directory of their own, since
`build.sh` builds without tests:

```bash
emcmake cmake -S . -B build-wasm-tests -DCMAKE_BUILD_TYPE=Release -DSSC_WASM_THREADS=ON -DSSC_BUILD_TESTS=ON
cmake --build build-wasm-tests --parallel
node build-wasm-tests/ssc_tests.js
```

All 88 test cases pass, the same as in the native build. Every method and
solve-for choice with its defaults gives the native answer to within 4 parts
in 10¹², the difference between two math libraries. CI runs these tests and
builds the page on every push.
