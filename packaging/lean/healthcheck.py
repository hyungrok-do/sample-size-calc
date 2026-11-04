import json
import os
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

request = Request('http://127.0.0.1:8080/api/lean/health',
                  headers={'Host': urlsplit(os.environ['SSC_LEAN_PUBLIC_ORIGIN']).netloc})
with urlopen(request, timeout=8) as response:
    if not json.load(response)['ready']:
        raise SystemExit(1)
