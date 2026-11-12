#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REMOTE_DIR="${REMOTE_DIR:-${ROOT_DIR}/../tucanos-remote-mdbook-test}"
TOOLCHAIN="${RUSTUP_TOOLCHAIN:-1.99.1}"
TMP_WORK_DIR="${TMPDIR:-${ROOT_DIR}/.tmp-rustdoc}"

if [[ ! -d "${REMOTE_DIR}/.git" ]]; then
  echo "Cloning remote tucanos repository into ${REMOTE_DIR}"
  git clone https://github.com/tucanos/tucanos "${REMOTE_DIR}"
fi

echo "Building tmesh/tucanos with Rust ${TOOLCHAIN}"
(
  cd "${REMOTE_DIR}"
  RUSTUP_TOOLCHAIN="${TOOLCHAIN}" cargo build -p tmesh -p tucanos
)

echo "Running mdbook test with linked remote crates"
mkdir -p "${TMP_WORK_DIR}"
(
  cd "${ROOT_DIR}"
  TMPDIR="${TMP_WORK_DIR}" \
  RUSTUP_TOOLCHAIN="${TOOLCHAIN}" \
  mdbook test . \
    -L "${REMOTE_DIR}/target/debug/deps" \
    -L "${REMOTE_DIR}/target/debug"
)

echo "mdbook test completed successfully"
