// What a static host such as GitHub Pages cannot do for the page, done in the browser:
//  - the two headers that make the page cross-origin isolated, so the engine can use every core (from the
//    visit after this worker is installed);
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
  'format.js?v=dev'];
// The only copies kept, each under the one address the page asks for it by: any other address (another query
// string on a visit, a link to a file with ?x= added) is answered from the network and not stored, so links
// cannot fill the site's storage with copies.
const PAGE = new URL('./', self.location).href;
const KEPT = new Set([...FILES, 'worker.js?engine=ssc_web_mt&v=dev', 'engine/ssc_web_mt.mjs?v=dev',
  'engine/ssc_web_mt.wasm?v=dev'].map((file) => new URL(file, self.location).href));

function keptAs(request) {
  const url = new URL(request.url);
  url.hash = '';
  if (request.mode !== 'navigate') return KEPT.has(url.href) ? url.href : null;
  // a visit to the page, whatever its query, is kept as the page; a visit to any other file is not
  const path = url.origin + url.pathname;
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

function isolated(response) {
  if (response.status === 0) return response;  // an opaque response cannot be changed
  const headers = new Headers(response.headers);
  headers.set('Cross-Origin-Opener-Policy', 'same-origin');
  headers.set('Cross-Origin-Embedder-Policy', 'require-corp');
  headers.set('Cross-Origin-Resource-Policy', 'same-origin');
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
      return isolated(response);
    } catch (error) {
      const cached = await caches.match(keptAs(request) ?? request, { ignoreSearch: true });
      if (!cached) throw error;
      return isolated(cached);
    }
  })());
});
