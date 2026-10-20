"""Serves the web page on this computer: python3 web/serve.py [port], then open http://localhost:8000.

The page is served cross-origin isolated (the two headers below), which lets the engine use every core; add
?engine=single to the address to try the one-thread engine that any web server can serve.
"""
import functools
import http.server
import sys
from pathlib import Path


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {
        **http.server.SimpleHTTPRequestHandler.extensions_map,
        ".js": "text/javascript",
        ".mjs": "text/javascript",
        ".wasm": "application/wasm",
    }

    def end_headers(self):
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()


def main() -> None:
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    here = Path(__file__).resolve().parent
    if not (here / "engine" / "ssc_web.mjs").is_file():
        sys.exit("The engine is not built yet: run web/build.sh first.")
    handler = functools.partial(Handler, directory=str(here))
    with http.server.ThreadingHTTPServer(("127.0.0.1", port), handler) as server:
        print(f"Serving the calculator at http://localhost:{port} (Ctrl+C to stop)")
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
