//! Host-side build helper for the plugin components.
//!
//! `wasm-tools component new` needs the MSVC linker to build on this machine, so the
//! same step is done here with the `wit-component` library it is built on.
//!
//! Usage:
//!   xtask component <core.wasm> <out.wasm>   wrap a core module (with wit-bindgen metadata) as a component
//!   xtask imports <core.wasm>                list a core module's imports (spot WASI / wasm-bindgen / env)

use anyhow::{Context, Result, bail};
use std::{env, fs};

fn main() -> Result<()> {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.iter().map(String::as_str).collect::<Vec<_>>().as_slice() {
        ["component", input, output] => component(input, output),
        ["imports", input] => imports(input),
        _ => bail!("usage: xtask component <core.wasm> <out.wasm> | xtask imports <core.wasm>"),
    }
}

fn component(input: &str, output: &str) -> Result<()> {
    let core = fs::read(input).with_context(|| format!("reading {input}"))?;
    let foreign = foreign_imports(&core)?;
    if !foreign.is_empty() {
        bail!("{input} imports modules the SolidCraft host does not provide:\n  {}", foreign.join("\n  "));
    }
    let bytes = wit_component::ComponentEncoder::default()
        .module(&core)?
        .validate(true)
        .encode()
        .context("encoding the component")?;
    fs::write(output, &bytes).with_context(|| format!("writing {output}"))?;
    println!("{output}: component, {} bytes (core {} bytes)", bytes.len(), core.len());
    Ok(())
}

fn imports(input: &str) -> Result<()> {
    let core = fs::read(input).with_context(|| format!("reading {input}"))?;
    for (module, name) in core_imports(&core)? {
        println!("{module} :: {name}");
    }
    Ok(())
}

/// Every (module, name) a core module imports.
fn core_imports(core: &[u8]) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for payload in wasmparser::Parser::new(0).parse_all(core) {
        if let wasmparser::Payload::ImportSection(section) = payload? {
            for import in section.into_imports() {
                let import = import?;
                out.push((import.module.to_string(), import.name.to_string()));
            }
        }
    }
    Ok(out)
}

/// Imports that are neither the plugin WIT (`craft:plugin/...`) nor the canonical-ABI shims.
fn foreign_imports(core: &[u8]) -> Result<Vec<String>> {
    Ok(core_imports(core)?
        .into_iter()
        .filter(|(module, _)| !module.starts_with("craft:plugin/") && !module.starts_with("[export]craft:plugin/") && module != "$root")
        .map(|(module, name)| format!("{module} :: {name}"))
        .collect())
}
