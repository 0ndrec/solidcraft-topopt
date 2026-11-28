//! The docked panel: its layout (read from the manifest, so `describe()` always matches what was
//! approved at install) and its data (the stored study).

use crate::craft::plugin::types::{ColorRole, Datum, PluginError, Value as Cell};
use crate::craft::plugin::ui::{self, *};
use serde_json::Value;
use topopt_study::study::fmt_n;
use topopt_study::solver::SolverView;
use topopt_study::{BodyInfo, Study};

const MANIFEST: &str = include_str!("../../topopt/plugin.json");
pub const PANEL: &str = "topopt";

fn s(v: &Value) -> String {
    v.as_str().unwrap_or_default().to_owned()
}
fn opt(v: &Value) -> Option<String> {
    v.as_str().map(str::to_owned)
}

fn button(v: &Value) -> ButtonSpec {
    ButtonSpec { label: s(&v["label"]), icon: opt(&v["icon"]), entry: s(&v["entry"]) }
}

fn widget(w: &Value) -> Widget {
    let (kind, v) = w.as_object().and_then(|o| o.iter().next()).map(|(k, v)| (k.as_str(), v)).unwrap_or(("text", &Value::Null));
    let ids = |v: &Value| v.as_array().map(|a| a.iter().filter_map(|x| x.as_u64()).map(|x| x as u32).collect()).unwrap_or_default();
    match kind {
        "group" => Widget::Group(ids(v)),
        "tabs" => Widget::Tabs(
            v.as_array().into_iter().flatten().map(|t| Tab { label: s(&t["label"]), child: t["child"].as_u64().unwrap_or(0) as u32 }).collect(),
        ),
        "split" => Widget::Split(SplitSpec {
            vertical: v["vertical"].as_bool().unwrap_or(false),
            first: v["first"].as_u64().unwrap_or(0) as u32,
            second: v["second"].as_u64().unwrap_or(0) as u32,
        }),
        "toolbar" => Widget::Toolbar(v.as_array().into_iter().flatten().map(button).collect()),
        "button" => Widget::Button(button(v)),
        "form" => Widget::Form(FormSpec {
            bind: s(&v["bind"]),
            fields: v["fields"]
                .as_array()
                .into_iter()
                .flatten()
                .map(|f| Field {
                    key: s(&f["key"]),
                    label: s(&f["label"]),
                    kind: match f["kind"].as_str() {
                        Some("length") => FieldKind::Length,
                        Some("angle") => FieldKind::Angle,
                        Some("text") => FieldKind::Text,
                        Some("toggle") => FieldKind::Toggle,
                        Some("choice") => FieldKind::Choice,
                        Some("selection") => FieldKind::Selection,
                        _ => FieldKind::Number,
                    },
                    choices: f["choices"].as_array().into_iter().flatten().map(s).collect(),
                    on_change: opt(&f["on_change"]),
                })
                .collect(),
        }),
        "table" => Widget::Table(TableSpec {
            bind: s(&v["bind"]),
            row_id: s(&v["row_id"]),
            columns: v["columns"]
                .as_array()
                .into_iter()
                .flatten()
                .map(|c| Column { key: s(&c["key"]), label: s(&c["label"]), unit: opt(&c["unit"]) })
                .collect(),
            select: match v["select"].as_str() {
                Some("single") => SelectionMode::Single,
                Some("multi") => SelectionMode::Multi,
                _ => SelectionMode::None,
            },
            on_select: opt(&v["on_select"]),
            on_edit: opt(&v["on_edit"]),
        }),
        "item_list" => Widget::ItemList(ListSpec { bind: s(&v["bind"]), on_select: opt(&v["on_select"]) }),
        "chart" => Widget::Chart(ChartSpec { kind: ChartKind::Line, series: vec![ChartSeries { label: String::new(), role: ColorRole::Series1, x: String::new(), y: String::new() }] }),
        "progress" => Widget::Progress(s(v)),
        "empty_state" => Widget::EmptyState(s(v)),
        _ => Widget::Text(s(v)),
    }
}

/// The panel declarations of `runtime.ui` in plugin.json.
pub fn declarations() -> Declarations {
    let m: Value = serde_json::from_str(MANIFEST).expect("plugin.json is valid JSON");
    let u = &m["runtime"]["ui"];
    let panels = u["panels"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|p| PanelSpec {
            id: s(&p["id"]),
            title: s(&p["title"]),
            zones: p["zones"]
                .as_array()
                .into_iter()
                .flatten()
                .map(|z| match z.as_str() {
                    Some("left") => DockZone::Left,
                    Some("bottom") => DockZone::Bottom,
                    _ => DockZone::Right,
                })
                .collect(),
            layout: Layout {
                root: p["layout"]["root"].as_u64().unwrap_or(0) as u32,
                nodes: p["layout"]["nodes"].as_array().into_iter().flatten().map(|n| Node { widget: widget(&n["widget"]) }).collect(),
            },
        })
        .collect();
    Declarations { ui_version: u["ui_version"].as_u64().unwrap_or(1) as u32, panels, workspaces: Vec::new() }
}

fn text(path: String, t: impl Into<String>) -> Datum {
    Datum { path, value: Cell::Text(t.into()) }
}
fn number(path: String, n: f64) -> Datum {
    Datum { path, value: Cell::Number(n) }
}
fn flag(path: String, b: bool) -> Datum {
    Datum { path, value: Cell::Flag(b) }
}

/// Rows of a two-column key/value table.
fn kv(bind: &str, rows: &[(&str, String)], out: &mut Vec<Datum>) {
    for (i, (k, v)) in rows.iter().enumerate() {
        out.push(text(format!("{bind}/{i}/key"), *k));
        out.push(text(format!("{bind}/{i}/value"), v.clone()));
    }
}

/// Replaces the panel's data with the study. `body`: the measured design body, when found.
pub fn show(study: &Study, body: Option<&BodyInfo>, solver: &SolverView) -> Result<(), PluginError> {
    let mut d = Vec::new();
    for (row, c) in study.conditions.iter().enumerate() {
        let (kind, value) = c.summary();
        d.push(text(format!("/items/{row}/id"), c.id().to_string()));
        d.push(text(format!("/items/{row}/kind"), kind));
        d.push(text(format!("/items/{row}/where"), c.where_()));
        d.push(text(format!("/items/{row}/value"), value));
    }
    let rows = topopt_study::status::rows_with(study, body, solver);
    // The first row ("Next") also stands above the tabs, so it is seen from every tab.
    if let Some((_, next)) = rows.first() {
        kv("/next", &[("Next", next.clone())], &mut d);
    }
    let refs: Vec<(&str, String)> = rows.iter().skip(1).map(|(k, v)| (k.as_str(), v.clone())).collect();
    kv("/status", &refs, &mut d);

    let st = &study.settings;
    for (k, v) in [
        ("volume_percent", (st.volfrac * 1000.0).round() / 10.0),
        ("element_size", st.element_size),
        ("max_elements", st.max_elements as f64),
        ("min_member", st.min_member),
        ("pad", st.pad),
        ("max_iterations", st.max_iterations as f64),
        ("time_limit", st.time_limit),
        ("penal", st.penal),
        ("beta_max", st.beta_max),
        ("move_limit", st.move_limit),
    ] {
        d.push(number(format!("/settings/{k}"), v));
    }
    d.push(flag("/settings/use_solver".into(), st.use_solver));
    d.push(flag("/settings/connect_supports".into(), st.connect_supports));
    let m = &study.material;
    d.push(text("/material/label".into(), m.label()));
    d.push(number("/material/e".into(), m.e));
    d.push(number("/material/nu".into(), m.nu));
    d.push(number("/material/density".into(), m.density));
    let ms = &study.mesh;
    d.push(number("/mesh/triangles".into(), ms.triangles as f64));
    d.push(number("/mesh/smoothing".into(), ms.smoothing as f64));
    d.push(number("/mesh/iso".into(), ms.iso));
    d.push(flag("/mesh/keep_faces".into(), ms.keep_faces));
    d.push(flag("/mesh/remove_islands".into(), ms.remove_islands));
    d.push(flag("/mesh/fill_cavities".into(), ms.fill_cavities));
    d.push(flag("/mesh/convert_solid".into(), ms.convert_solid));
    // One toggle per smoother, ticked for the current one (see the manifest generator).
    for (id, _) in topopt_study::SMOOTHERS {
        d.push(flag(format!("/mesh/smoother_{id}"), ms.smoother == id));
    }
    d.push(flag("/mesh/preserve_volume".into(), ms.preserve_volume));
    d.push(flag("/mesh/sharp_edges".into(), ms.sharp_edges));
    d.push(number("/mesh/upsample".into(), ms.upsample as f64));
    d.push(number("/mesh/trim_percent".into(), ms.trim_percent));
    d.push(flag("/mesh/fillets".into(), ms.fillets));
    d.push(number("/mesh/fillet_radius".into(), ms.fillet_radius));

    let mut rows: Vec<(&str, String)> = Vec::new();
    let ready = solver.ready.as_ref().and_then(|f| f.result.clone());
    match ready.as_ref().or(study.result.as_ref()) {
        Some(r) => {
            rows.push(("Mass", format!("{:.1} g", r.mass_g)));
            if r.max_displacement_mm > 0.0 {
                rows.push(("Max deflection", format!("{:.3} mm", r.max_displacement_mm)));
                rows.push(("Stress (von Mises)", format!("99 % below {:.1} MPa, peak {:.1} MPa", r.stress_p99_mpa, r.stress_max_mpa)));
                if r.cases.len() > 1 {
                    for c in &r.cases {
                        rows.push(("Load case", format!("{}: {:.3} mm, {:.1} MPa (99 %)", c.case, c.max_displacement_mm, c.stress_p99_mpa)));
                    }
                }
            }
            rows.push(("Volume", format!("{:.1} %, {} mm³", r.volume_fraction * 100.0, fmt_n(r.mesh_volume_mm3.round()))));
            rows.push(("Compliance", format!("{:.4} N·mm", r.compliance)));
            rows.push(("Run", format!("{}, {} it., {}", r.outcome, r.iterations, topopt_study::status::fmt_duration(r.seconds))));
            if !r.solver.is_empty() {
                rows.push(("Solved on", r.solver.clone()));
            }
            rows.push(("Cells", format!("{} of {:.2} mm", fmt_n(r.elements as f64), r.h)));
            rows.push(("Mesh", format!("{} triangles", fmt_n(r.triangles as f64))));
            rows.push(("Holes through", r.genus.to_string()));
            if r.trimmed > 0.0 {
                rows.push(("Trimmed", format!("{:.1} % of the material", r.trimmed * 100.0)));
            }
            rows.push(("File", format!("results/{}", r.stl)));
            rows.push(("In the design", if r.inserted { "yes".into() } else { "no: click Create body".into() }));
            for w in &r.warnings {
                rows.push(("Warning", w.clone()));
            }
        }
        None => rows.push(("—", "No result yet".into())),
    }
    kv("/result", &rows, &mut d);

    ui::clear_data(PANEL)?;
    ui::set_data(PANEL, &d)?;
    Ok(())
}
