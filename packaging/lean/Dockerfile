# syntax=docker/dockerfile:1
# Build the immutable compiler worker and the complete three-tab web application.
FROM python:3.12-slim-bookworm AS lean-worker
ARG TARGETARCH
RUN test "$TARGETARCH" = amd64
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential ca-certificates curl git zstd libgmp-dev \
    && rm -rf /var/lib/apt/lists/*
COPY lean/lean-toolchain /tmp/lean-toolchain
RUN version=$(sed 's|leanprover/lean4:v||' /tmp/lean-toolchain) \
    && mkdir -p /opt/lean \
    && curl --fail --location --retry 3 \
       "https://github.com/leanprover/lean4/releases/download/v${version}/lean-${version}-linux.tar.zst" \
       | tar --zstd -x -C /opt/lean --strip-components=1
ENV PATH="/opt/lean/bin:${PATH}" PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
WORKDIR /opt/compiler/lean
COPY lean/lakefile.toml lean/lake-manifest.json lean/lean-toolchain ./
RUN lake exe cache get Mathlib.Probability.Distributions.Gaussian Mathlib.Probability.CDF \
    Mathlib.Probability.Variance Mathlib.Tactic
COPY lean/ ./
COPY scripts/check_lean_audit.py /opt/compiler/scripts/
RUN lake build DesignCompiler && lake env lean Audit.lean > /tmp/audit.txt \
    && python /opt/compiler/scripts/check_lean_audit.py /tmp/audit.txt \
    && rm /tmp/audit.txt
COPY scripts/lean_equations.py scripts/lean_cpp.py scripts/lean_study_cpp.py /opt/compiler/scripts/
COPY packaging/lean/freeze_environment.py /opt/compiler/freeze_environment.py
RUN lake env python /opt/compiler/freeze_environment.py /opt/compiler/lean-env.json
COPY packaging/lean/worker.py /opt/compiler/worker.py
LABEL ssc.lean.worker="1"
USER 65532:65532
ENTRYPOINT ["python", "/opt/compiler/worker.py"]

FROM emscripten/emsdk:6.0.13 AS web-build
WORKDIR /src
COPY CMakeLists.txt LICENSE THIRD_PARTY_NOTICES.md ./
COPY core/ core/
COPY web/build.sh web/build.sh
RUN web/build.sh
COPY lean/ lean/
COPY packaging/lean-wasm/ packaging/lean-wasm/
COPY scripts/build_lean_wasm.py scripts/build_lean_wasm.py
RUN python3 scripts/build_lean_wasm.py --output /browser-compiler
COPY web/ web/
RUN cp -R /browser-compiler web/lean-wasm
ARG SITE_VERSION=container
RUN web/assemble.sh /site "$SITE_VERSION"

FROM docker:28-cli AS docker-client
FROM python:3.12-slim-bookworm AS gateway
COPY --from=docker-client /usr/local/bin/docker /usr/local/bin/docker
RUN pip install --no-cache-dir gunicorn==26.2.0 && docker --version
WORKDIR /app
COPY packaging/lean/app.py packaging/lean/containers.py packaging/lean/healthcheck.py ./
COPY --from=web-build /site/ /srv/site/
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
EXPOSE 8080
HEALTHCHECK --interval=30s --timeout=10s --start-period=15s CMD ["python", "healthcheck.py"]
# A single process owns the capacity/cancellation registry. Threads keep health
# and cancellation responsive while at most SSC_LEAN_MAX_JOBS workers compile.
CMD ["gunicorn", "app:create_app()", "--bind=0.0.0.0:8080", "--workers=1", "--threads=16", "--timeout=120", "--graceful-timeout=90", "--no-control-socket", "--access-logfile=-", "--error-logfile=-", "--limit-request-line=2048", "--limit-request-fields=40", "--limit-request-field_size=4096"]
