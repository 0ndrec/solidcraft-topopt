"""Writes plugins/topopt/plugin.json: entries (dialogs), panel layout and permissions of the
TopOpt plugin. The plugin's ui-provider.describe() reads the same file, so they cannot drift.

    python -I scripts/gen_topopt_icons.py && python -I scripts/gen_topopt_manifest.py

Every entry shows in SolidCraft's Plugins group and menu, with its menu-slot name as its title and
its icon. So there are few entries, each a readable action: the panel's buttons run the dialog
entries themselves (which then open their form), and one `panel` entry takes the panel's field
changes and row double-clicks (and, run from the menu, opens the panel).
"""
import json
import pathlib

MATERIALS = [
    ("aluminium-6061", "Aluminium 6061-T6"), ("aluminium-7075", "Aluminium 7075-T6"), ("steel", "Steel (structural)"),
    ("stainless-316", "Stainless steel 316"), ("titanium", "Titanium Ti-6Al-4V"), ("pla", "PLA (3D print)"),
    ("petg", "PETG (3D print)"), ("abs", "ABS (3D print)"), ("nylon", "Nylon PA12"),
]


def sel(title, accept, desc, min_items=None, max_items=None):
    p = {"type": "array", "title": title, "description": desc, "x-solidcraft-selection": {"version": 1, "accept": [accept]}}
    if min_items is not None:
        p["minItems"] = min_items
    if max_items is not None:
        p["maxItems"] = max_items
    return p


def num(title, default, lo, hi, desc, unit=None, integer=False):
    p = {"type": "integer" if integer else "number", "title": title, "default": default, "minimum": lo, "maximum": hi, "description": desc}
    if unit:
        p["unit"] = unit
    return p


def flag(title, default, desc):
    return {"type": "boolean", "title": title, "default": default, "description": desc}


def named(title, icon, schema):
    """A readable title (the entry's menu-slot name: menus, ribbon and dialog show it) and an icon."""
    schema = dict(schema)
    schema["x-solidcraft-icon"] = f"icons/{icon}.svg"
    return {"slots": [{"kind": "menu", "name": title}], "input_schema": schema}


EMPTY = {"type": "object", "properties": {}}
PANEL_EVENT = {"type": "object", "additionalProperties": True, "properties": {}}

entries = [
    {"name": "design-body", "effect": "edit", **named("Design body", "design_body", {
        "type": "object", "required": ["body"], "properties": {
            "body": sel("Design body", "body", "The solid that bounds the part: material is only placed inside it.", 1, 1),
            "keep_out": sel("Keep empty", "body", "Optional: bodies whose space must stay free (bolt heads, cables, other parts)."),
            "material": {"type": "string", "title": "Material", "default": "aluminium-6061",
                         "oneOf": [{"const": i, "title": t} for i, t in MATERIALS] + [{"const": "custom", "title": "Custom (set E, nu, density in the panel)"}]},
            "volume_percent": num("Keep", 30.0, 2.0, 95.0, "How much of the design body's volume the result may use, in percent.", "%"),
            "element_size": num("Element size", 0.0, 0.0, 10000.0, "Size of one voxel. 0 = from the element budget.", "mm"),
            "max_elements": num("Element budget", 40000, 1000, 400000, "More elements give finer detail and take longer (40 000: a few minutes).", integer=True),
            "min_member": num("Thinnest member", 0.0, 0.0, 10000.0, "Smallest strut or wall thickness. 0 = about three elements.", "mm"),
        }})},
    {"name": "add-support", "effect": "edit", **named("Add support", "add_support", {
        "type": "object", "required": ["faces"], "properties": {
            "faces": sel("Faces", "face", "Faces of the design body that are bolted, clamped or resting on something.", 1),
            "fix_x": flag("Hold in X", True, "Untick to let the faces slide along X."),
            "fix_y": flag("Hold in Y", True, "Untick to let the faces slide along Y."),
            "fix_z": flag("Hold in Z", True, "Untick to let the faces slide along Z."),
            "replace": num("Item to change", 0, 0, 100000, "Filled in when you double-click a row; 0 adds a new support.", integer=True),
            "delete": flag("Delete this item", False, "With an item number: remove it instead of changing it."),
        }})},
    {"name": "add-load", "effect": "edit", **named("Add load", "add_load", {
        "type": "object", "required": ["faces"], "properties": {
            "faces": sel("Faces", "face", "Faces of the design body the force acts on; it is spread evenly over them.", 1),
            "mode": {"type": "string", "title": "Direction", "default": "vector", "oneOf": [
                {"const": "vector", "title": "Along X, Y, Z (enter Fx, Fy, Fz)"},
                {"const": "normal", "title": "Pushing into the faces (enter Magnitude)"}]},
            "fx": num("Fx", 0.0, -1e9, 1e9, "Total force along X.", "N"),
            "fy": num("Fy", 0.0, -1e9, 1e9, "Total force along Y.", "N"),
            "fz": num("Fz", -1000.0, -1e9, 1e9, "Total force along Z (negative = down).", "N"),
            "magnitude": num("Magnitude", 1000.0, 0.0, 1e9, "For Pushing into the faces.", "N"),
            "case": num("Load case", 1, 1, 8, "Loads with the same number act together; different numbers are separate situations the part must survive.", integer=True),
            "replace": num("Item to change", 0, 0, 100000, "Filled in when you double-click a row; 0 adds a new load.", integer=True),
            "delete": flag("Delete this item", False, "With an item number: remove it instead of changing it."),
        }})},
    {"name": "keep-face", "effect": "edit", **named("Keep surface", "keep_face", {
        "type": "object", "required": ["faces"], "properties": {
            "faces": sel("Faces", "face", "Faces of the design body that must stay exactly as they are: mating, sealing or cosmetic faces.", 1),
            "thickness": num("Solid depth", 0.0, 0.0, 10000.0, "Material kept under the faces. 0 = the pad thickness (about 1.5 elements).", "mm"),
            "replace": num("Item to change", 0, 0, 100000, "Filled in when you double-click a row; 0 adds a new kept surface.", integer=True),
            "delete": flag("Delete this item", False, "With an item number: remove it instead of changing it."),
        }})},
    {"name": "optimize", "job": True, **named("Optimize", "optimize", EMPTY)},
    {"name": "create-body", "effect": "edit", **named("Create body", "create_body", EMPTY)},
    {"name": "remesh", "job": True, **named("Re-mesh", "remesh", EMPTY)},
    {"name": "delete-item", "effect": "edit", **named("Delete item", "delete_item", {
        "type": "object", "required": ["id"], "properties": {"id": num("Item", 1, 1, 100000, "The number in the # column of the panel.", integer=True)}})},
    {"name": "clear-all", "effect": "edit", **named("Clear supports and loads", "clear_all", EMPTY)},
    {"name": "refresh", "effect": "read", **named("Refresh TopOpt panel", "refresh", EMPTY)},
    # Opens the panel; also takes the panel's field changes and row double-clicks.
    {"name": "panel", "effect": "edit", **named("TopOpt panel", "panel", PANEL_EVENT)},
]


def btn(label, entry):
    return {"label": label, "icon": None, "entry": entry}


def field(key, label, kind, choices=None):
    return {"key": key, "label": label, "kind": kind, "choices": choices or [], "on_change": "panel"}


def col(key, label, unit=None):
    return {"key": key, "label": label, "unit": unit}


def table(bind, columns, select="none", on_edit=None):
    return {"table": {"bind": bind, "row_id": columns[0]["key"], "columns": columns, "select": select, "on_select": None, "on_edit": on_edit}}


nodes = []


def node(widget):
    nodes.append({"widget": widget})
    return len(nodes) - 1


# A Next line, then tabs: every page fits the panel (the host does not scroll panels; tables
# scroll themselves).
root = node(None)  # filled below
next_row = node(table("/next", [col("key", "Step"), col("value", "Next")]))
intro = node({"text": "Pick the design body, then the faces held (supports) and pushed (loads); optionally faces that must not change."})
setup_tools = node({"toolbar": [btn("Design body…", "design-body"), btn("Add support…", "add-support"), btn("Add load…", "add-load"),
                                btn("Keep surface…", "keep-face")]})
items = node(table("/items", [col("id", "#"), col("kind", "Type"), col("where", "Where"), col("value", "Value")], "single", "panel"))
hint = node({"text": "Double-click a row to change or delete it."})
item_tools = node({"toolbar": [btn("Delete…", "delete-item"), btn("Clear all", "clear-all")]})
setup = node({"group": [intro, setup_tools, items, hint, item_tools]})

run_tools = node({"toolbar": [btn("Optimize", "optimize"), btn("Create body", "create-body"), btn("Re-mesh", "remesh"), btn("Refresh", "refresh")]})
run_note = node({"text": "Optimize runs in the background (operations panel, bottom right). With TopOpt Solver running it uses every CPU core."})
status = node(table("/status", [col("key", "Study"), col("value", "Value")]))
result = node(table("/result", [col("key", "Result"), col("value", "Value")]))
run = node({"group": [run_tools, result, run_note, status]})

basic = node({"form": {"bind": "/settings", "fields": [
    field("volume_percent", "Keep, % of the design body", "number"),
    field("element_size", "Element size (0 = from budget)", "length"),
    field("max_elements", "Element budget", "number"),
    field("min_member", "Thinnest member (0 = auto)", "length"),
    field("connect_supports", "Join every support to the part", "toggle"),
    field("max_iterations", "Max iterations", "number"),
    field("time_limit", "Time limit, s", "number")]}})
material = node({"form": {"bind": "/material", "fields": [
    field("label", "Material", "choice", [t for _, t in MATERIALS] + ["Custom"]),
    field("e", "Young's modulus E, MPa", "number"),
    field("nu", "Poisson ratio", "number"),
    field("density", "Density, g/cm3", "number")]}})
mesh = node({"form": {"bind": "/mesh", "fields": [
    field("triangles", "Triangles (up to 50 000)", "number"),
    field("smoothing", "Smoothing passes", "number"),
    # The smoother as three toggles that act as radio buttons: SolidCraft's panel dropdowns
    # (0.49-0.50) never send a change event, so a "choice" field cannot be changed.
    field("smoother_bilateral", "Smoother: bilateral (keeps creases)", "toggle"),
    field("smoother_taubin", "Smoother: Taubin (fastest)", "toggle"),
    field("smoother_hc", "Smoother: HC (keeps volume)", "toggle"),
    field("upsample", "Finer surface (1 = off, 2 = round members)", "number"),
    field("keep_faces", "Keep the body's faces exact", "toggle"),
    field("sharp_edges", "Keep its sharp edges sharp", "toggle"),
    field("fillets", "Round the joints to supports and loads", "toggle"),
    field("remove_islands", "Remove floating pieces", "toggle"),
    field("fill_cavities", "Fill closed cavities", "toggle"),
    field("convert_solid", "Create body makes a solid", "toggle")]}})
adv_settings = node({"form": {"bind": "/settings", "fields": [
    field("use_solver", "Use TopOpt Solver when it is running", "toggle"),
    field("pad", "Solid pad under supports/loads/kept faces (0 = auto)", "length"),
    field("penal", "SIMP penalty", "number"),
    field("beta_max", "Final sharpness (beta)", "number"),
    field("move_limit", "Move limit per iteration", "number")]}})
adv_mesh = node({"form": {"bind": "/mesh", "fields": [field("iso", "Surface density threshold", "number"),
                                                    field("preserve_volume", "Keep the volume while smoothing", "toggle"),
                                                    field("fillet_radius", "Fillet radius (0 = 1.5 elements)", "length"),
                                                    field("trim_percent", "Trim dead weight costing up to % stiffness (0 = off)", "number")]}})
advanced = node({"group": [adv_settings, adv_mesh]})
tabs = node({"tabs": [{"label": "Setup", "child": setup}, {"label": "Run", "child": run}, {"label": "Settings", "child": basic},
                      {"label": "Material", "child": material}, {"label": "Mesh", "child": mesh}, {"label": "Advanced", "child": advanced}]})
nodes[root] = {"widget": {"group": [next_row, tabs]}}

manifest = {
    "id": "topopt",
    "name": "TopOpt",
    "version": "0.12.0",
    "description": "Topology optimization: pick a design body, supports and loads; get a lightweight, stiff part as a mesh or solid.",
    "capabilities": ["commands@0.2", "ui@0.1", "storage@0.1", "forms@0.1"],
    "runtime": {
        "module": "topopt.wasm",
        "data_format": 1,
        "commands": ["storage.get", "storage.set", "inspect.measure", "model.faces", "view.pick", "model.mesh", "file.insert_mesh", "mesh.convert"],
        "entries": entries,
        "ui": {"ui_version": 1, "panels": [{"id": "topopt", "title": "TopOpt", "zones": ["right", "left", "bottom"],
                                            "layout": {"root": root, "nodes": nodes}}], "workspaces": []},
        "permissions": {"command_groups": ["file"], "files": [{"name": "results", "read": True, "write": True}], "clock": True},
        "limits": {"fuel": 10_000_000_000_000, "memory_mib": 1024, "deadline_ms": 30000},
    },
}

out = pathlib.Path(__file__).resolve().parent.parent / "plugins" / "topopt" / "plugin.json"
# Compact: the manifest limit is 16 KiB.
text = json.dumps(manifest, separators=(",", ":"), ensure_ascii=False) + "\n"
assert len(text.encode()) < 16 * 1024, f"manifest is {len(text.encode())} bytes, over 16 KiB"
out.write_text(text, encoding="utf-8")
print(out, len(text.encode()), "bytes,", len(entries), "entries")
