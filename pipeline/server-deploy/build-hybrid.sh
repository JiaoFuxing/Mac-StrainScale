#!/bin/bash
# Build hybrid on the server (optional; the prebuilt static hybrid usually works as-is).
# Option A (default): glibc dynamic build, needs a Rust toolchain on the server.
# Option B (musl static): needs cargo-zigbuild + zig (see README).
set -e
cd "$(dirname "$0")/hybrid-rust"
if [ "${1:-}" = "musl" ]; then
    cargo-zigbuild build --release --target x86_64-unknown-linux-musl
    echo "OK: target/x86_64-unknown-linux-musl/release/hybrid (static)"
else
    cargo build --release
    echo "OK: target/release/hybrid (glibc dynamic)"
fi
