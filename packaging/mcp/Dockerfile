# The MCP server, hosted: Streamable HTTP on $PORT (default 8000) at /mcp, /healthz for health checks.
#   docker build -f packaging/mcp/Dockerfile -t ssc-mcp .
#   docker run -p 8000:8000 -e SSC_MCP_ALLOWED_HOSTS=mcp.example.org ssc-mcp
# Only the engine and the MCP SDK are installed (no Qt): the server needs nothing else.
FROM python:3.12-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends build-essential cmake ninja-build \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
RUN pip wheel --no-deps -w /wheels .

FROM python:3.12-slim
RUN useradd --create-home ssc
COPY --from=build /wheels /wheels
RUN pip install --no-cache-dir --no-deps /wheels/*.whl && pip install --no-cache-dir "mcp>=2,<3" && rm -rf /wheels
USER ssc
ENV PORT=8000
EXPOSE 8000
HEALTHCHECK CMD python -c "import urllib.request, os; urllib.request.urlopen(f'http://127.0.0.1:{os.environ[\"PORT\"]}/healthz')"
CMD ["sh", "-c", "exec ssc-mcp --transport http --host 0.0.0.0 --port ${PORT}"]
