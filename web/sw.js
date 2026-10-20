// What a static host such as GitHub Pages cannot do for the page, done in the browser:
//  - the two headers that make the page cross-origin isolated, so the engine can use every core (from the
//    visit after this worker is installed);
//  - a copy of each of the page's files, used when the network is not there, so the page and its engine
//    start offline.
// The network comes first, so a new version is used as soon as it is published.
const CACHE = 'ssc-1';
// fetched at installation: the first visit loaded them before this worker could keep a copy (the multi-core
// engine is kept when it is first used, from the next visit)
const FILES = ['./', 'app.js', 'engine.js', 'worker.js', 'engine/ssc_web.mjs', 'engine/ssc_web.wasm'];

self.addEventListener('install', (event) => {
  self.skipWaiting();
  event.waitUntil(caches.open(CACHE).then((cache) => cache.addAll(FILES)).catch(() => {}));
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
      const response = await fetch(request);
      if (response.ok) {
        const copy = response.clone();
        event.waitUntil(caches.open(CACHE).then((cache) => cache.put(request, copy)));
      }
      return isolated(response);
    } catch (error) {
      const cached = await caches.match(request, { ignoreSearch: true });
      if (!cached) throw error;
      return isolated(cached);
    }
  })());
});
