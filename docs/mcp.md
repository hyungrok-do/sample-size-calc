# The calculator for language models (MCP server)

`ssc-mcp` makes the calculator available to language-model clients (Claude,
ChatGPT, Codex, coding agents) through the
[Model Context Protocol](https://modelcontextprotocol.io). The model chooses
the method and fills in the inputs; the engine does the arithmetic, the same
engine and the same validated results as the web page and the command line.

The server is an outer layer: like the command line, it calls the engine
through the Python bindings (`ssc.core`). Nothing else goes through the
server, and nothing else needs the MCP package.

```
C++ engine  ←  ssc.core  ←  ssc-calc · ssc-mcp
C++ engine (WebAssembly)  ←  the web page
```

## What the model gets

Three tools, the loop a person follows on the web page or the command line:

| Tool | What it does |
|---|---|
| `search_methods` | Methods that fit a description of a study ("two independent means", "logrank non-inferiority", "cluster randomized proportions"), best first; an empty description lists all 71 |
| `describe_method` | A method's solve options and inputs: ids, choices, defaults, allowed values, and when each input applies; also an example and the method's caveats |
| `calculate` | Sample size, power or the smallest detectable effect, for one or many scenarios, with the plain-language statement of each, the dropout-inflated sizes, the notes, and warnings about anything in the request that was not used |

All three are read-only: they compute from their arguments and touch nothing
else (`readOnlyHint: true`, `destructiveHint: false`, `idempotentHint: true`,
`openWorldHint: false`). Results are compact JSON.

What makes it safe for a model to use:

- **Nothing is dropped silently.** An input that the chosen options leave out
  (`n2` with equal groups), a choice not available in this mode (`alloc=sep`
  when solving for the sample size) or an unknown input name is reported,
  with when the input applies, so the model can correct the request.
- **Errors say what to change**: an effect on the wrong side of a one-sided
  test, a target power at or below alpha, an unknown choice with the valid
  ones listed.
- **Simulation methods** return their quick analytic approximation unless
  `simulate=true` asks for the Monte Carlo simulation, which takes longer.

## Install

The server needs the engine and the MCP SDK (Python 3.10 or later):

```bash
pip install "sample-size-calculator[mcp] @ git+https://github.com/hyungrok-do/sample-size-calc.git"
# or, from a clone:  pip install ".[mcp]"
```

This builds the C++ engine, so it needs CMake and a C++17 compiler (see
[Building and testing from source](../README.md#building-and-testing-from-source)).

## Use it from a local client

`ssc-mcp` talks over standard input and output, which local clients start
themselves.

**Claude Code**

```bash
claude mcp add sample-size -- ssc-mcp
```

**Claude Desktop**: in `claude_desktop_config.json` (Settings → Developer →
Edit Config):

```json
{
  "mcpServers": {
    "sample-size": { "command": "ssc-mcp" }
  }
}
```

**Codex**: in `~/.codex/config.toml`:

```toml
[mcp_servers.sample-size]
command = "ssc-mcp"
```

Other clients take the same command. If `ssc-mcp` is not on the client's
`PATH`, give its full path (`which ssc-mcp`).

## Host it (Streamable HTTP)

ChatGPT, claude.ai and the directories below connect to a server on the
internet. `ssc-mcp --transport http` serves the MCP endpoint at `/mcp`, a
health check at `/healthz`, and the ChatGPT domain-verification token (below)
at `/.well-known/openai-apps-challenge`.

```bash
docker build -f packaging/mcp/Dockerfile -t ssc-mcp .
docker run -p 8000:8000 -e SSC_MCP_ALLOWED_HOSTS=mcp.example.org ssc-mcp
```

The image holds only the engine and the MCP SDK. Put it behind HTTPS (any
container host, such as Cloud Run, Fly.io or Render, or a reverse proxy with a
certificate). It needs no sign-in: the server computes from the numbers it is
sent and stores nothing.

| Variable | Default (hosted) | What it sets |
|---|---|---|
| `PORT` | 8000 | the port |
| `SSC_MCP_ALLOWED_HOSTS` | — | the public host name(s), comma-separated: requests for other `Host` headers are refused (DNS-rebinding protection) |
| `SSC_MCP_TIME_LIMIT` | 60 | seconds per calculation before it is stopped |
| `SSC_MCP_MAX_SCENARIOS` | 200 | input combinations per request (at most the engine's 20,000) |
| `SSC_MCP_MAX_REPS` | 5000 | Monte Carlo replicates per scenario |
| `SSC_MCP_THREADS` | 2 | cores per calculation |
| `SSC_OPENAI_APPS_CHALLENGE` | — | the token ChatGPT asks to find at `/.well-known/openai-apps-challenge` |

Over stdio the limits are generous (600 seconds, 20,000 scenarios, every
core), since the server runs on the user's own computer.

## Get listed in the Claude and ChatGPT directories

Listing is done by the project's owner with their own accounts; this
repository has what the reviews check. Both directories take a **hosted**
server (previous section): an `https://` URL ending in `/mcp`.

### Claude: Connectors Directory

Submit at [claude.ai/directory/manage](https://claude.ai/directory/manage)
→ **Submit new** → **MCP connector** (any paid Claude plan). Requirements,
from [Anthropic's submission guide](https://claude.com/docs/connectors/building/submission):

| Requirement | Status |
|---|---|
| Remote server over HTTPS | host the image above |
| Every tool has a `title` and `readOnlyHint` or `destructiveHint` | done (all three tools) |
| Authentication: OAuth for a user's account, or none for public data | none: the server holds no user data |
| Tested as a custom connector in Claude | add the URL under Settings → Connectors → Add custom connector and run each tool |
| Documentation URL | this page |
| Privacy policy URL | [docs/privacy.md](privacy.md) |
| Support contact, icon | the issues page; `packaging/mcp/icon.png` (512 × 512) |
| Reviewer instructions | the test prompts below (no account is needed) |

Anthropic scans the submission and lists it as a Community connector; some
submissions also get a person's review. Local servers are no longer accepted
as desktop extensions; a local server reaches the directory only inside a
plugin.

### ChatGPT (and Codex): Plugins directory

OpenAI lists MCP servers as plugins, shared by ChatGPT and Codex, submitted
from the OpenAI Platform by a verified organization. See
[Remote MCP server review requirements](https://developers.openai.com/plugins/deploy/app-review)
and the [plugin guidelines](https://developers.openai.com/apps-sdk/app-submission-guidelines).

| Requirement | Status |
|---|---|
| Production HTTPS MCP endpoint (no tunnels or staging URLs) | host the image above |
| Domain verification at `/.well-known/openai-apps-challenge` | set `SSC_OPENAI_APPS_CHALLENGE` to the token OpenAI gives |
| `readOnlyHint`, `openWorldHint` and `destructiveHint` set on every tool, with a justification for each | done; justification: the tools compute from their arguments, read and change nothing, and reach nothing outside the server |
| Five positive and three negative test cases | below |
| Name, descriptions, icon, screenshots, privacy policy, terms, support contact | privacy policy above; the rest at submission |

To try it before submitting, add the URL as a custom connector in ChatGPT's
developer mode.

### The official MCP Registry

[registry.modelcontextprotocol.io](https://registry.modelcontextprotocol.io)
lists servers for every client. `packaging/mcp/server.json` describes this one
for a PyPI release; the README carries the `mcp-name` line the registry checks
on PyPI. Once the package is on PyPI (or, for the hosted server, after adding
a `remotes` entry with its URL), publish with the registry's `mcp-publisher`
tool:

```bash
mcp-publisher login github      # proves the io.github.hyungrok-do/ namespace
mcp-publisher publish packaging/mcp/server.json
```

Check the file with `mcp-publisher validate` first, as the schema evolves.

## Test prompts for reviewers

Each prompt should lead the model to `search_methods`, then `describe_method`,
then `calculate`, and to answer with the statement.

**Positive**

1. "How many patients per group do I need to detect a 5 mmHg difference in
   blood pressure, with a standard deviation of 12 and 90% power?" → two-sample
   t-test: 123 per group, 246 in total (90.2% power).
2. "What power does a trial with 150 patients per arm have to detect response
   rates of 40% versus 25%?" → two proportions: 79.6% power.
3. "Hazard ratio 0.7, control median survival 12 months, 12 months of accrual
   and 12 of follow-up: how many patients for 80% power?" → logrank: 215 per
   group, about 248 events.
4. "A cluster-randomized trial with 20 patients per clinic, ICC 0.05, comparing
   30% versus 20%: how many clinics per arm for 80% power?" → 29 clinics per
   arm, 1,160 patients.
5. "With 50 per group, what is the smallest standardized difference a t-test
   detects with 80% power?" → 0.566 standard deviations.

**Negative**

1. "Find the sample size for 4% power." → refused: the target power must
   exceed alpha.
2. "One-sided test that the treatment raises the mean, expecting it to lower
   the mean by 0.5." → refused: the difference is on the side opposite to the
   one-sided alternative.
3. "What is the weather in Paris?" → the tools are not used.
