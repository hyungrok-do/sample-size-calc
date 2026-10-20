# Packaging

`mcp/` holds what the MCP server (`ssc-mcp`) needs to be hosted and listed:

| File | Purpose |
|---|---|
| `mcp/Dockerfile` | An image that serves `ssc-mcp` over Streamable HTTP: `docker build -f packaging/mcp/Dockerfile -t ssc-mcp .` |
| `mcp/server.json` | The server's entry in the MCP Registry (`mcp-publisher publish packaging/mcp/server.json`) |
| `mcp/icon.png` | The icon for directory listings (512 × 512) |

[`docs/mcp.md`](../docs/mcp.md) covers the clients, hosting and listing. The
Python package itself (the engine, `ssc-calc` and `ssc-mcp`) is built by
`pip` from `pyproject.toml`; the web page is built by
[`web/build.sh`](../web/build.sh).
