#!/usr/bin/env bash
# Cross-compile axia CLI for the targets this host can produce and pack them.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SOLVER="$ROOT/solver"
DIST="$ROOT/dist"
VER="$(sed -n 's/^version = "\(.*\)"/\1/p' "$SOLVER/Cargo.toml" | head -1)"
TRIPLE_HOST="$(rustc -vV | awk '/host:/{print $2}')"

mkdir -p "$DIST"
echo "Axia FEM ${VER}  host=${TRIPLE_HOST}"

pack() {
  local target="$1"
  local bin="$2"
  local src="$SOLVER/target/${target}/release/${bin}"
  if [[ ! -f "$src" ]]; then
    echo "skip pack ${target}: ${src} missing" >&2
    return 1
  fi
  local stage="axia-${VER}-${target}"
  rm -rf "$DIST/$stage"
  mkdir -p "$DIST/$stage/examples"
  cp "$src" "$DIST/$stage/"
  cp "$ROOT/README.md" "$ROOT/LICENSE" "$DIST/$stage/"
  cp "$ROOT/examples/"*.inp "$DIST/$stage/examples/" 2>/dev/null || true
  if [[ "$bin" == *.exe ]]; then
    python3 - <<PY
import zipfile, os
dist = "$DIST"
stage = "$stage"
zpath = os.path.join(dist, stage + ".zip")
with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(os.path.join(dist, stage)):
        for f in files:
            p = os.path.join(root, f)
            z.write(p, os.path.relpath(p, dist))
print("packed ", zpath)
PY
  else
    tar -C "$DIST" -czf "$DIST/${stage}.tar.gz" "$stage"
    echo "packed  $DIST/${stage}.tar.gz"
  fi
}

build_one() {
  local target="$1"
  echo
  echo "==> $target"
  local cmd=(cargo build --release --bin axia --features cli --target "$target")
  if command -v cargo-zigbuild >/dev/null 2>&1 && command -v zig >/dev/null 2>&1; then
    cmd=(cargo zigbuild --release --bin axia --features cli --target "$target")
  fi
  echo "    ${cmd[*]}"
  if (cd "$SOLVER" && "${cmd[@]}"); then
    if [[ "$target" == *windows* ]]; then
      pack "$target" axia.exe
    else
      pack "$target" axia
    fi
  else
    echo "BUILD FAILED: $target" >&2
    return 1
  fi
}

TARGETS=("$@")
if [[ ${#TARGETS[@]} -eq 0 ]]; then
  TARGETS=("$TRIPLE_HOST")
  rustup target list --installed | grep -qx 'x86_64-unknown-linux-musl' && TARGETS+=("x86_64-unknown-linux-musl")
  rustup target list --installed | grep -qx 'aarch64-unknown-linux-gnu' && TARGETS+=("aarch64-unknown-linux-gnu")
  rustup target list --installed | grep -qx 'x86_64-pc-windows-gnu' && TARGETS+=("x86_64-pc-windows-gnu")
  rustup target list --installed | grep -qx 'aarch64-apple-darwin' && TARGETS+=("aarch64-apple-darwin")
  rustup target list --installed | grep -qx 'x86_64-apple-darwin' && TARGETS+=("x86_64-apple-darwin")
fi

ok=()
fail=()
for t in "${TARGETS[@]}"; do
  if build_one "$t"; then
    ok+=("$t")
  else
    fail+=("$t")
  fi
done

echo
echo "---- checksums ----"
(cd "$DIST" && sha256sum axia-${VER}-*.tar.gz axia-${VER}-*.zip 2>/dev/null || true) | tee "$DIST/SHA256SUMS"

echo
echo "built: ${ok[*]}"
if [[ ${#fail[@]} -gt 0 ]]; then
  echo "failed: ${fail[*]}" >&2
  exit 1
fi
