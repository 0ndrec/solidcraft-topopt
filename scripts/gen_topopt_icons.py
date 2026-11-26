"""Writes the TopOpt plugin's icons to plugins/topopt/icons/*.svg, in SolidCraft's icon subset:
viewBox 0 0 64 64, absolute M L C Q A Z paths, palette colours or currentColor, round strokes.
The style follows the app's model icons: an isometric block with light / mid / dark faces and
ink outlines, blue for what the action creates, orange for what acts on the part.

    python -I scripts/gen_topopt_icons.py
"""
import pathlib
import re

INK, P1, P2, P3 = "#34465b", "#e2e8ef", "#b8c3cf", "#8e9faf"
A1, A2, A3 = "#b9ddf8", "#328bd3", "#236ead"
O1, O2, O3 = "#f1d2af", "#c06f00", "#a45000"
W = "2.55"


def pts(*p):
    return "M" + "L".join(f"{x:g} {y:g}" for x, y in p) + "Z"


def stroke(d, color=INK, width=W):
    return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>'


def fill(d, color):
    return f'<path d="{d}" fill="{color}"/>'


def block(cx, top, w, h, depth, faces=(P1, P2, P3)):
    """An isometric block: top face centred at cx, `top` its highest y; w half-width, h height."""
    t, r, b, l = (cx, top), (cx + w, top + depth), (cx, top + 2 * depth), (cx - w, top + depth)
    l2, b2, r2 = (l[0], l[1] + h), (b[0], b[1] + h), (r[0], r[1] + h)
    parts = [fill(pts(l, b, b2, l2), faces[1]), fill(pts(b, r, r2, b2), faces[2]), fill(pts(t, r, b, l), faces[0])]
    parts.append(stroke(pts(t, r, r2, b2, l2, l)))
    parts.append(stroke(f"M{l[0]:g} {l[1]:g}L{b[0]:g} {b[1]:g}L{r[0]:g} {r[1]:g}"))
    parts.append(stroke(f"M{b[0]:g} {b[1]:g}L{b2[0]:g} {b2[1]:g}"))
    return parts


ICONS = {
    # The design space: a blue block.
    "design_body": ("Design body", block(32, 6, 25, 26, 13, (A1, A2, A3))),
    # A block resting on a hatched ground.
    "add_support": ("Add support", block(32, 4, 22, 18, 11, (P1, P2, P3)) + [
        stroke("M8 52L56 52", O3, "3.4"),
        stroke("M14 52L9 60M24 52L19 60M34 52L29 60M44 52L39 60M54 52L49 60", O2, "2.6"),
        fill(pts((22, 44), (32, 49), (42, 44), (32, 52)), O2),
    ]),
    # A force arrow pushing down on a block.
    "add_load": ("Add load", block(32, 26, 24, 16, 12, (P1, P2, P3)) + [
        stroke("M32 4L32 26", O2, "4.2"),
        fill(pts((23, 22), (41, 22), (32, 34)), O2),
    ]),
    # A block whose top face is kept: the face in blue, a padlock on it.
    "keep_face": ("Keep surface", block(28, 8, 23, 20, 12, (A2, P2, P3)) + [
        stroke("M43 44L43 39A6 6 0 0 1 55 39L55 44", O3, "3.2"),
        fill("M39 44L59 44L59 58L39 58Z", O2),
        fill("M47.5 48L50.5 48L50.5 54L47.5 54Z", "#ffffff"),
    ]),
    # The optimized part: a blue truss over the faint outline of its block.
    "optimize": ("Optimize", [
        stroke(pts((32, 6), (57, 19), (57, 45), (32, 58), (7, 45), (7, 19)), "#6f879c", "2"),
        stroke("M9 44L32 14L55 44L32 56Z", A3, "4.4"),
        stroke("M32 14L32 56M20 29L32 42L44 29", A2, "3.2"),
    ]),
    # A new body: a blue block with a plus.
    "create_body": ("Create body", block(28, 10, 22, 22, 11, (A1, A2, A3)) + [
        '<circle cx="49" cy="49" r="12" fill="#c06f00"/>',
        stroke("M49 42L49 56M42 49L56 49", "#ffffff", "3.4"),
    ]),
    # A triangle mesh.
    "remesh": ("Re-mesh", [
        fill(pts((32, 6), (58, 54), (6, 54)), A1),
        stroke(pts((32, 6), (58, 54), (6, 54)), A3, "3"),
        stroke("M19 30L45 30L32 54Z", A3, "2.4"),
    ]),
    # A list item with a cross.
    "delete_item": ("Delete item", [
        stroke("M8 16L34 16M8 30L28 30M8 44L22 44", "currentColor", "4"),
        '<circle cx="45" cy="42" r="13" fill="#c06f00"/>',
        stroke("M39.5 36.5L50.5 47.5M50.5 36.5L39.5 47.5", "#ffffff", "3.4"),
    ]),
    # A list wiped clear.
    "clear_all": ("Clear all", [
        stroke("M8 14L56 14M8 28L56 28M8 42L36 42", "currentColor", "4"),
        stroke("M42 38L58 54M58 38L42 54", O2, "4.4"),
    ]),
    # Circular arrows.
    "refresh": ("Refresh", [
        stroke("M50 24A20 20 0 0 0 14 26", "currentColor", "4.25"),
        fill(pts((53, 12), (54, 30), (37, 26)), "currentColor"),
        stroke("M14 40A20 20 0 0 0 50 38", "currentColor", "4.25"),
        fill(pts((11, 52), (10, 34), (27, 38)), "currentColor"),
    ]),
    # A docked panel.
    "panel": ("TopOpt panel", [
        stroke("M10 9L54 9Q59 9 59 14L59 50Q59 55 54 55L10 55Q5 55 5 50L5 14Q5 9 10 9Z", "currentColor", "4.25"),
        fill("M36 13L54 13Q55 13 55 14L55 50Q55 51 54 51L36 51Z", A2),
        stroke("M40 21L50 21M40 29L50 29M40 37L47 37", "#ffffff", "2.6"),
    ]),
}

PATH_RE = re.compile(r"[MLCQAZ0-9 .\-]+")
out = pathlib.Path(__file__).resolve().parent.parent / "plugins" / "topopt" / "icons"
out.mkdir(parents=True, exist_ok=True)
for name, (title, parts) in ICONS.items():
    body = "\n  ".join(parts)
    svg = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">\n  <title>{title}</title>\n  {body}\n</svg>\n'
    for d in re.findall(r' d="([^"]*)"', svg):
        assert PATH_RE.fullmatch(d), (name, d)
    (out / f"{name}.svg").write_text(svg, encoding="utf-8")
print(f"{len(ICONS)} icons in {out}")
