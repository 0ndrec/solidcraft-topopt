#!/usr/bin/env bash
# Build TopOpt Solver (native, release) and put it with its launcher into build/solver/.
#
#   scripts/build-solver.sh
#
# A running solver keeps its exe locked on Windows: it is renamed out of the way first.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cargo +stable-x86_64-pc-windows-gnu build --release -p topopt-solver --manifest-path "$root/Cargo.toml"
mkdir -p "$root/build/solver" "$root/build/results"
exe="$root/build/solver/topopt-solver.exe"
if [ -f "$exe" ] && ! cp "$root/target/release/topopt-solver.exe" "$exe" 2>/dev/null; then
    mv "$exe" "$exe.old" && cp "$root/target/release/topopt-solver.exe" "$exe"
fi
[ -f "$exe" ] || cp "$root/target/release/topopt-solver.exe" "$exe"
cp "$root/tools/topopt-solver/start-topopt-solver.cmd" "$root/build/solver/"
"$exe" version
