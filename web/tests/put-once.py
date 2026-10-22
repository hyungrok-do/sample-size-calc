import sys
p = sys.argv[1] + "/sw.js"
s = open(p).read()
old = """      const key = keptAs(request);
      if (response.ok && key) {
        const copy = response.clone();
        event.waitUntil(caches.open(CACHE).then((cache) => cache.put(key, copy)));
      }"""
new = """      const key = keptAs(request);
      // a file is the same under its versioned address for good: kept once; the page itself every time
      if (response.ok && key && (key === PAGE || !(await caches.match(key)))) {
        const copy = response.clone();
        event.waitUntil(caches.open(CACHE).then((cache) => cache.put(key, copy)));
      }"""
assert old in s, "sw.js changed"
open(p, "w").write(s.replace(old, new))
