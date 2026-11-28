"""Renders a mesh (PLY or STL) to a PNG with a transparent background, for the documentation.
A small software rasterizer (numpy, z-buffer, smooth shading, 2× supersampling): no GPU and no
plotting library needed.

    tools/reference/.venv/Scripts/python scripts/render_mesh.py MESH OUT.png [--size W H]
        [--yaw DEG] [--pitch DEG] [--color RRGGBB] [--up z|y] [--compare MESH2]...

With --compare (repeatable), each further mesh is rendered with the same camera and framing and
put to the right of MESH (a before/after picture, or a row of variants).
"""
import argparse
import struct
import zlib

import numpy as np
import trimesh


def write_png(path, rgba):
    h, w, _ = rgba.shape
    raw = b"".join(b"\x00" + rgba[y].tobytes() for y in range(h))

    def chunk(kind, data):
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)

    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0)))
        f.write(chunk(b"IDAT", zlib.compress(raw, 9)))
        f.write(chunk(b"IEND", b""))


def rotation(yaw, pitch):
    """World → camera: turn about Z by yaw, then tilt by pitch. Camera looks along −z'."""
    a, b = np.radians(yaw), np.radians(pitch)
    rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    # Z up in the world becomes y up on screen; then tilt down by pitch.
    to_screen = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
    rx = np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)], [0, np.sin(b), np.cos(b)]])
    return rx @ to_screen @ rz


def render(mesh, width, height, yaw, pitch, color, frame=None):
    """`frame`: the vertices that set the framing (default: the mesh's own)."""
    # Normals are not averaged across sharp edges: a large flat face whose corners it shares with
    # a steep wall would otherwise be shaded in blotches.
    mesh = mesh.smooth_shaded
    ss = 2
    W, H = width * ss, height * ss
    v = np.asarray(mesh.vertices, dtype=np.float64)
    R = rotation(yaw, pitch)
    p = v @ R.T
    f = (np.asarray(frame, dtype=np.float64) if frame is not None else v) @ R.T
    centre = (f.min(0) + f.max(0)) / 2  # centre the projection, not the world box
    p -= centre
    f -= centre
    n = np.asarray(mesh.vertex_normals) @ R.T
    span = max((f[:, 0].max() - f[:, 0].min()) / W, (f[:, 1].max() - f[:, 1].min()) / H) * 1.08
    x = p[:, 0] / span + W / 2
    y = H / 2 - p[:, 1] / span
    z = p[:, 2]
    faces = np.asarray(mesh.faces)
    # Back faces are skipped: on a closed mesh they are hidden anyway, and on thin plates their
    # depth is so close to the front faces' that they would speckle through.
    fn = np.cross(p[faces[:, 1]] - p[faces[:, 0]], p[faces[:, 2]] - p[faces[:, 0]])
    front = faces[fn[:, 2] > 0]

    depth = np.full((H, W), -np.inf)
    normal = np.zeros((H, W, 3))
    for t in front:
        xs, ys = x[t], y[t]
        x0, x1 = int(max(np.floor(xs.min()), 0)), int(min(np.ceil(xs.max()), W - 1))
        y0, y1 = int(max(np.floor(ys.min()), 0)), int(min(np.ceil(ys.max()), H - 1))
        if x0 > x1 or y0 > y1:
            continue
        px, py = np.meshgrid(np.arange(x0, x1 + 1) + 0.5, np.arange(y0, y1 + 1) + 0.5)
        d = (ys[1] - ys[2]) * (xs[0] - xs[2]) + (xs[2] - xs[1]) * (ys[0] - ys[2])
        if abs(d) < 1e-12:
            continue
        l0 = ((ys[1] - ys[2]) * (px - xs[2]) + (xs[2] - xs[1]) * (py - ys[2])) / d
        l1 = ((ys[2] - ys[0]) * (px - xs[2]) + (xs[0] - xs[2]) * (py - ys[2])) / d
        l2 = 1 - l0 - l1
        inside = (l0 >= 0) & (l1 >= 0) & (l2 >= 0)
        if not inside.any():
            continue
        zz = l0 * z[t[0]] + l1 * z[t[1]] + l2 * z[t[2]]
        region = depth[y0:y1 + 1, x0:x1 + 1]
        win = inside & (zz > region)
        if not win.any():
            continue
        region[win] = zz[win]
        nn = l0[..., None] * n[t[0]] + l1[..., None] * n[t[1]] + l2[..., None] * n[t[2]]
        normal[y0:y1 + 1, x0:x1 + 1][win] = nn[win]

    covered = np.isfinite(depth)
    nrm = normal / np.maximum(np.linalg.norm(normal, axis=2, keepdims=True), 1e-12)
    # Camera looks along −z: a visible face's normal has z > 0.
    lights = [(np.array([-0.4, 0.6, 0.7]), 0.75), (np.array([0.6, 0.2, 0.5]), 0.30), (np.array([0.0, -0.8, 0.3]), 0.12)]
    shade = np.full((H, W), 0.22)
    for d, k in lights:
        d = d / np.linalg.norm(d)
        shade += k * np.clip((nrm * d).sum(2), 0, None)
    view = np.array([0, 0, 1.0])
    half = np.array([-0.4, 0.6, 1.7])
    half /= np.linalg.norm(half)
    spec = np.clip((nrm * half).sum(2), 0, None) ** 40 * 0.35
    rim = (1 - np.clip((nrm * view).sum(2), 0, 1)) ** 3 * 0.15
    base = np.array([int(color[i:i + 2], 16) / 255 for i in (0, 2, 4)])
    rgb = np.clip(base[None, None] * shade[..., None] + spec[..., None] + rim[..., None], 0, 1)

    # Outline where depth jumps (silhouettes and creases seen edge-on).
    dz = np.zeros((H, W))
    finite = np.where(covered, depth, depth[covered].min() - (np.ptp(depth[covered]) + 1))
    dz[1:, :] = np.maximum(dz[1:, :], np.abs(np.diff(finite, axis=0)))
    dz[:, 1:] = np.maximum(dz[:, 1:], np.abs(np.diff(finite, axis=1)))
    edge = dz > 4 * span
    rgb[edge] = rgb[edge] * 0.45
    alpha = covered.astype(float)
    alpha[edge & ~covered] = 0

    rgba = np.concatenate([rgb * alpha[..., None], alpha[..., None]], axis=2)
    rgba = rgba.reshape(height, ss, width, ss, 4).mean(axis=(1, 3))
    a = rgba[..., 3:4]
    rgba[..., :3] = np.where(a > 0, rgba[..., :3] / np.maximum(a, 1e-6), 0)
    return (np.clip(rgba, 0, 1) * 255 + 0.5).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("out")
    ap.add_argument("--size", nargs=2, type=int, default=[900, 600])
    ap.add_argument("--yaw", type=float, default=-35)
    ap.add_argument("--pitch", type=float, default=28)
    ap.add_argument("--color", default="8fc1ea")
    ap.add_argument("--up", choices=["z", "y"], default="z")
    ap.add_argument("--compare", action="append", default=[], help="another mesh, rendered alike to the right (repeat for more)")
    a = ap.parse_args()

    def load(path):
        mesh = trimesh.load(path, force="mesh", process=True)
        if a.up == "y":
            mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [1, 0, 0]))
        # Split the vertices at sharp creases, so smooth shading does not blur them.
        return trimesh.graph.smooth_shade(mesh, angle=np.radians(60))

    mesh = load(a.mesh)
    if a.compare:
        meshes = [mesh] + [load(c) for c in a.compare]
        frame = np.vstack([m.vertices for m in meshes])
        gap = np.zeros((a.size[1], a.size[0] // 20, 4), dtype=np.uint8)
        tiles = []
        for m in meshes:
            if tiles:
                tiles.append(gap)
            tiles.append(render(m, a.size[0], a.size[1], a.yaw, a.pitch, a.color, frame))
        write_png(a.out, np.concatenate(tiles, axis=1))
    else:
        write_png(a.out, render(mesh, a.size[0], a.size[1], a.yaw, a.pitch, a.color))
    print(a.out, len(mesh.faces), "triangles")


if __name__ == "__main__":
    main()
