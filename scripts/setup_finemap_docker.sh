#!/usr/bin/env bash
set -euo pipefail

# Build a local FINEMAP wrapper from official FINEMAP binary tarball.
# This is designed for Apple Silicon/macOS where native FINEMAP binary linkage
# can be problematic.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS_DIR="$ROOT_DIR/tools/finemap"
DOCKER_DIR="$TOOLS_DIR/docker"
BIN_DIR="$ROOT_DIR/tools/bin"

MAC_TAR="$TOOLS_DIR/finemap_v1.4.2_MacOSX.tgz"
UNIX_TAR="$TOOLS_DIR/finemap_v1.4.2_x86_64.tgz"
MAC_URL="https://www.christianbenner.com/finemap_v1.4.2_MacOSX.tgz"
UNIX_URL="https://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz"

mkdir -p "$TOOLS_DIR" "$DOCKER_DIR" "$BIN_DIR"

echo ">>> command -v docker"
command -v docker >/dev/null

if [[ ! -f "$MAC_TAR" ]]; then
  echo ">>> curl -k -fL -o $MAC_TAR $MAC_URL"
  curl -k -fL -o "$MAC_TAR" "$MAC_URL"
fi

if [[ ! -f "$UNIX_TAR" ]]; then
  echo ">>> curl -k -fL -o $UNIX_TAR $UNIX_URL"
  curl -k -fL -o "$UNIX_TAR" "$UNIX_URL"
fi

echo ">>> tar -tzf $UNIX_TAR | head -n 10"
tar -tzf "$UNIX_TAR" | head -n 10

cat > "$DOCKER_DIR/Dockerfile" <<'EOF'
FROM --platform=linux/amd64 ubuntu:22.04
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates zstd libgomp1 && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/finemap
COPY finemap_v1.4.2_x86_64.tgz /opt/finemap/
RUN tar -xzf finemap_v1.4.2_x86_64.tgz && chmod +x /opt/finemap/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64
ENTRYPOINT ["/opt/finemap/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64"]
EOF

cp "$UNIX_TAR" "$DOCKER_DIR/"

echo ">>> docker build --platform linux/amd64 -t ckd-finemap:1.4.2 $DOCKER_DIR"
docker build --platform linux/amd64 -t ckd-finemap:1.4.2 "$DOCKER_DIR"

cat > "$BIN_DIR/finemap" <<'EOF'
#!/bin/zsh
set -euo pipefail
IMAGE="ckd-finemap:1.4.2"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORKDIR="${PWD}"
if [[ "$WORKDIR" != "$PROJECT_ROOT"* ]]; then
  WORKDIR="$PROJECT_ROOT"
fi
exec docker run --rm --platform linux/amd64 \
  -v "$PROJECT_ROOT":"$PROJECT_ROOT" \
  -w "$WORKDIR" \
  "$IMAGE" "$@"
EOF
chmod +x "$BIN_DIR/finemap"

echo ">>> PATH=\"$BIN_DIR:\$PATH\" finemap --help | head -n 20"
PATH="$BIN_DIR:$PATH" finemap --help | head -n 20

echo "FINEMAP setup complete."
echo "Add this to your shell for benchmarks:"
echo "  export PATH=\"$BIN_DIR:\$PATH\""
