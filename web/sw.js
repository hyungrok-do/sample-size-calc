// What a static host such as GitHub Pages cannot do for the page, done in the browser:
//  - the two headers that make the page cross-origin isolated, so the engine can use every core (from the
//    visit after this worker is installed);
//  - the page's Content-Security-Policy as a header, which also reaches the engine's worker and its threads (the
//    <meta> in the page covers the page alone) and forbids other sites to frame the page (a <meta> cannot);
//  - a copy of each of the page's files, used when the network is not there, so the page and its engine
//    start offline.
// The network comes first, so a new version is used as soon as it is published.  Every file but the page is
// asked for with the version in its address (?v=, stamped by web/stamp.sh when the site is published), so a
// browser never mixes a new page with an older file it kept.
const CACHE = 'ssc-dev';
// fetched at installation: the first visit loaded them before this worker could keep a copy (the multi-core
// engine is kept when it is first used, from the next visit)
const FILES = ['./', 'app.js?v=dev', 'engine.js?v=dev', 'worker.js?engine=ssc_web&v=dev', 'engine/ssc_web.mjs?v=dev',
  'engine/ssc_web.wasm?v=dev', 'guide.js?v=dev', 'guide-data.js?v=dev', 'guide-state.js?v=dev', 'selection.js?v=dev',
  'format.js?v=dev', 'dag.html?v=dev', 'dag.css?v=dev', 'dag-entry.js?v=dev', 'dag.js?v=dev',
  'dag-model.js?v=dev', 'dag-view.js?v=dev', 'estimand-model.js?v=dev',
  'power-settings.js?v=dev', 'power-model.js?v=dev', 'dag-power.js?v=dev'];
// The only copies kept, each under the one address the page asks for it by: any other address (another query
// string on a visit, a link to a file with ?x= added) is answered from the network and not stored, so links
// cannot fill the site's storage with copies.
const PAGE = new URL('./', self.location).href;
const DAG_PAGE = new URL('dag.html?v=dev', self.location).href;
const KEPT = new Set([...FILES, 'worker.js?engine=ssc_web_mt&v=dev', 'engine/ssc_web_mt.mjs?v=dev',
  'engine/ssc_web_mt.wasm?v=dev'].map((file) => new URL(file, self.location).href));

function keptAs(request) {
  const url = new URL(request.url);
  url.hash = '';
  if (request.mode !== 'navigate') return KEPT.has(url.href) ? url.href : null;
  // a visit to the page, whatever its query, is kept as the page; a visit to any other file is not
  const path = url.origin + url.pathname;
  if (path === DAG_PAGE.split('?')[0]) return DAG_PAGE;
  return path === PAGE || path === `${PAGE}index.html` ? PAGE : null;
}

self.addEventListener('install', (event) => {
  self.skipWaiting();
  const fresh = FILES.map((file) => new Request(file, { cache: 'reload' }));  // not an older copy the browser kept
  event.waitUntil(caches.open(CACHE).then((cache) => cache.addAll(fresh)).catch(() => {}));
});

self.addEventListener('activate', (event) => {
  event.waitUntil((async () => {
    for (const key of await caches.keys()) if (key !== CACHE) await caches.delete(key);
    await self.clients.claim();
  })());
});

// The policy of the <meta> in index.html, plus frame-ancestors (npm test checks that the two agree).
const POLICY = "default-src 'none'; script-src 'self' 'wasm-unsafe-eval' 'unsafe-eval'; "
  + "style-src 'sha256-oUMCQLhHROcrKQ/eem2baHepydt9L9W1ITo5OyK16ro='; connect-src 'self'; worker-src 'self'; "
  + "img-src 'self'; manifest-src 'self'; base-uri 'none'; form-action 'none'; frame-ancestors 'none'";

// The independent DAG page uses its own stylesheet and the shared WebAssembly engine.
const DAG_POLICY = "default-src 'none'; script-src 'self' 'wasm-unsafe-eval'; style-src 'self'; connect-src 'self'; worker-src 'self'; "
  + "img-src 'self'; base-uri 'none'; form-action 'none'; frame-ancestors 'none'";

function isolated(response, request) {
  if (response.status === 0) return response;  // an opaque response cannot be changed
  const headers = new Headers(response.headers);
  headers.set('Cross-Origin-Opener-Policy', 'same-origin');
  headers.set('Cross-Origin-Embedder-Policy', 'require-corp');
  headers.set('Cross-Origin-Resource-Policy', 'same-origin');
  headers.set('Content-Security-Policy', POLICY);
  if (new URL(request.url).pathname === new URL(DAG_PAGE).pathname) headers.set('Content-Security-Policy', DAG_POLICY);
  return new Response(response.body, { status: response.status, statusText: response.statusText, headers });
}

self.addEventListener('fetch', (event) => {
  const request = event.request;
  if (request.method !== 'GET' || new URL(request.url).origin !== self.location.origin) return;
  event.respondWith((async () => {
    try {
      // the page itself is checked with the server every time: it names the version of every other file
      const response = await (request.mode === 'navigate' ? fetch(request.url, { cache: 'no-cache' }) : fetch(request));
      const key = keptAs(request);
      if (response.ok && key) {
        const copy = response.clone();
        event.waitUntil(caches.open(CACHE).then((cache) => cache.put(key, copy)));
      }
      return isolated(response, request);
    } catch (error) {
      const cached = await caches.match(keptAs(request) ?? request, { ignoreSearch: true });
      if (!cached) throw error;
      return isolated(cached, request);
    }
  })());
});
