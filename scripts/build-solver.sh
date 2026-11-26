#!/usr/bin/env bash
# Build TopOpt Solver and TopOpt CUDA Solver (native, release) and put them with their launchers
# into build/solver/.
#
#   scripts/build-solver.sh
#
# The CUDA solver needs no CUDA toolkit to build: its GPU kernels are checked in as PTX
# (scripts/build-ptx.sh rebuilds them). Without an NVIDIA GPU it solves on the CPU.
# A running solver keeps its exe locked on Windows: it is renamed out of the way first.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cargo +stable-x86_64-pc-windows-gnu build --release -p topopt-solver -p topopt-cuda-solver --manifest-path "$root/Cargo.toml"
mkdir -p "$root/build/solver" "$root/build/results"
for name in topopt-solver topopt-cuda-solver; do
    exe="$root/build/solver/$name.exe"
    if [ -f "$exe" ] && ! cp "$root/target/release/$name.exe" "$exe" 2>/dev/null; then
        mv "$exe" "$exe.old" && cp "$root/target/release/$name.exe" "$exe"
    fi
    [ -f "$exe" ] || cp "$root/target/release/$name.exe" "$exe"
    "$exe" version
done
cp "$root/tools/topopt-solver/start-topopt-solver.cmd" "$root/tools/topopt-cuda-solver/start-topopt-cuda-solver.cmd" "$root/build/solver/"
