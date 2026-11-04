# Deploy the Lean sample-size compiler

This deployment serves the complete calculator, DAG simulator, and Lean tab at
one HTTPS origin. The browser sends Lean source to `/api/lean/compile`; a fresh
worker compiles it with Lean 4.19.0 and the repository's pinned Mathlib, executes
the generated C++, and returns the derivation and numerical result. Users need
only a browser. The GitHub Pages site is still static and cannot run this backend.
Module and native-library paths are frozen when the image is built. Requests run
the installed Lean executable directly, with no Lake/Git initialization or
dependency downloads at runtime.

## Start on a Linux host

Use an x86-64 Linux host with Docker Engine and Docker Compose v2, a DNS name
pointing to the host, and inbound ports 80 and 443. Budget about 8 GB RAM and
40 GB free disk for building and running the default two-worker configuration,
including temporary Docker build layers and image-export copies.
The compiler and web assets are built from this checkout, without a prepublished
image or cloud-provider dependency.

```sh
cp packaging/lean/.env.example packaging/lean/.env
docker run --rm -it caddy:2.10.2-alpine caddy hash-password
```

Edit `packaging/lean/.env`: set `SSC_DOMAIN`, set your username, and paste the
generated password hash as the single-quoted `SSC_PASSWORD_HASH` value. Do not
put the plaintext password in that file. Set `SSC_VERSION` to the source commit
or another unique release label.

From the repository root:

```sh
docker compose --env-file packaging/lean/.env -f packaging/lean/compose.yaml build worker-image compiler
docker compose --env-file packaging/lean/.env -f packaging/lean/compose.yaml up -d compiler https
```

Open `https://YOUR_DOMAIN/lean.html`, sign in with the configured username and
password, and select **Compile & calculate**. Caddy provisions and renews TLS
certificates. The default repeated-measures study returns **51 participants per
arm, 102 total**. The header says **hosted compilation** when the backend is ready.
The other two tabs include their compiled single-threaded and multithreaded WASM
engines in the same deployment.

## Runtime boundaries

- The gateway serves static files and coordinates jobs. It does not execute Lean
  source or generated C++ in its own process/container.
- Each job gets a fresh worker with no network, no host bind mounts, a read-only
  root filesystem, UID 65532, dropped capabilities, and a private temporary
  filesystem. Neither the Docker socket nor login credentials enter workers.
- Defaults are two concurrent jobs, two CPUs and 2 GB memory per job, 128 processes,
  256 MB temporary storage, and a 60-second compiler deadline. The gateway also
  enforces a wall-clock deadline including worker startup and destroys workers on
  completion, cancellation, or failure. It removes leftover workers on restart.
- The gateway has Docker control privileges. Run it on a dedicated compiler host,
  keep its port private, and expose only the authenticated HTTPS proxy. The Compose
  file does this by publishing ports only on Caddy. A container is not a separate
  kernel: for broader untrusted access, configure a sandboxed OCI runtime such as
  gVisor with `SSC_RUNTIME`, or an equivalent VM-backed deployment.
- Source and results are transient. The editor retains source in the user's
  browser, and the normal download controls save modules, C++, and result JSON.
  Authentication identifies access to this deployment; it is not a multi-tenant
  account or project-storage system.

The worker uses the same compiler as local development. Container isolation does
not extend the mathematical verification claims: arithmetic/theorem checks and
floating-point/asymptotic power retain the distinctions described in
[the compiler documentation](../../docs/lean-compiler.md).

## Configuration and operation

| Variable | Default | Purpose |
| --- | --- | --- |
| `SSC_DOMAIN` | Required | Public DNS name for HTTPS |
| `SSC_USERNAME` | `researcher` | Login name |
| `SSC_PASSWORD_HASH` | Required | Caddy bcrypt password hash |
| `SSC_VERSION` | `local` | Image and static-asset release label |
| `SSC_MAX_JOBS` | `2` | Concurrent workers, 1–8 |
| `SSC_JOB_TIMEOUT` | `60` | Compiler seconds, 1–300 |
| `SSC_JOB_MEMORY_MB` | `2048` | Worker memory, 512–16384 MB |
| `SSC_JOB_CPUS` | `2` | Worker CPU quota and Lean threads, 1–16 |
| `SSC_DEPLOYMENT` | `sample-size-lean` | Worker ownership label, unique per controller |
| `SSC_RUNTIME` | Docker default | Optional installed OCI runtime name |

Run exactly one gateway process per deployment name. The gateway has a single
in-memory job/cancellation registry with multiple HTTP threads. Excess compilation
requests receive HTTP 429 with `Retry-After`, while health and cancellation remain
available. To scale across hosts, place a job queue and shared ownership registry
in front of separately isolated workers rather than increasing Gunicorn processes.

The authenticated `GET /api/lean/health` endpoint returns readiness and
`"execution":"hosted"`. `POST /api/lean/compile` takes
`{"source":"...", "jobId":"a-client-uuid"}`. `DELETE /api/lean/compile/JOB_ID`
cancels that job. Requests must come from the configured origin. The static site
and API are deployed together, so no CORS relaxation or browser API key is needed.

Inspect logs with:

```sh
docker compose --env-file packaging/lean/.env -f packaging/lean/compose.yaml logs --tail=100 compiler https
```

To update, pull the desired source revision, choose a new `SSC_VERSION`, and repeat
the build/start commands. The previous tagged images remain available for rollback.
Stop with the same Compose command followed by `down`. Keep the Caddy volumes to
retain certificates; add `--volumes` only when intentionally removing them.

## Verification

The `Lean deployment` workflow builds both real images, tests namespace/resource
isolation, compiles statistical studies in fresh workers, checks cancellation and
timeout cleanup, starts the HTTPS/authentication stack, and runs the browser suite
against that deployed stack. It checks both bundled WASM engines as well.

```sh
python scripts/check_lean_deployment.py
# After building the worker image:
python scripts/check_lean_deployment.py --docker
```

The smaller `--native` check tests the one-request worker protocol using an
already-built local Lean toolchain. It does not claim to test container isolation.

References: [Docker runtime constraints](https://docs.docker.com/engine/containers/run/),
[no-network workers](https://docs.docker.com/engine/network/drivers/none/), and
[Caddy authentication](https://caddyserver.com/docs/caddyfile/directives/basic_auth).
