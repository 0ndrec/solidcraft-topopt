#!/usr/bin/env bash
# Build one plugin package: plugins/<id>-src (Rust crate) -> plugins/<id>/<id>.wasm (component),
# then `solidcraft-cli plugin check`, and pack plugins/<id> (without tests/) into build/<id>.zip.
#
#   scripts/build-plugin.sh topopt
#
# Uses the self-contained GNU toolchain: this machine has no MSVC linker, which host-side build
# scripts and proc-macros would need under the default msvc toolchain.
set -euo pipefail
id="${1:?plugin id}"
root="$(cd "$(dirname "$0")/.." && pwd)"
tc="+stable-x86_64-pc-windows-gnu"
crate="$(sed -n 's/^name = "\(.*\)"/\1/p' "$root/plugins/$id-src/Cargo.toml" | head -1)"

cargo $tc build --release --target wasm32-unknown-unknown -p "$crate" --manifest-path "$root/Cargo.toml"
(cd "$root/xtask" && cargo $tc build --release -q)
core="$root/target/wasm32-unknown-unknown/release/${crate//-/_}.wasm"
"$root/xtask/target/release/xtask.exe" component "$core" "$root/plugins/$id/$id.wasm"
solidcraft-cli plugin check "$root/plugins/$id"

mkdir -p "$root/build"
rm -f "$root/build/$id.zip"
python -I - "$root/plugins/$id" "$root/build/$id.zip" <<'PY'
import sys, zipfile, pathlib
src, out = pathlib.Path(sys.argv[1]), sys.argv[2]
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for p in sorted(src.rglob("*")):
        rel = p.relative_to(src)
        if p.is_file() and rel.parts[0] != "tests" and not any(part.startswith(".") for part in rel.parts):
            z.write(p, rel.as_posix())
print(out)
PY
