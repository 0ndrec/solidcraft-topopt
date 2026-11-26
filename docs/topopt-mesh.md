# topopt-mesh — from a density field to a closed mesh

`crates/topopt-mesh` turns a density field (values 0…1 on a regular lattice) into a triangle mesh that is guaranteed **closed, 2-manifold, outward-oriented and free of self-intersections**: SolidCraft accepts it as a closed mesh body and converts it to a solid; trimesh, slicers and Blender see it as watertight.

The crate knows nothing about optimization: it needs only the field, so it can surface any scalar field (an SDF, another optimizer's result, a CT scan).

## The pipeline: `extract`

```rust
pub fn extract(field: Field, anchors: Option<&[bool]>, surface: Option<&Projector>, opts: &Options)
    -> Result<(Mesh, Report), String>
```

1. **Field** (`field`): optional blur; the border zeroed so the surface closes; values too close to `iso` moved off it by `nudge`; then **component clean-up** with the exact connectivity of the contour (14 neighbours, matching the tetrahedral split): islands without an anchor (`anchors`, usually the pads under supports and loads; without anchors, all but the largest) are removed and closed cavities filled.
2. **Contour** (`mt`): marching tetrahedra on the Freudenthal split (6 tetrahedra per cube, consistent diagonals between neighbours). Inside a tetrahedron linear interpolation gives a flat triangle or a convex quad, so the result is watertight, manifold and self-intersection-free **by construction**, without the ambiguous cases of marching cubes. One vertex per edge (hashed by its two samples); orientation from a table built from exact corner coordinates.
3. **Smoothing** (`smooth`): Taubin λ|μ, without shrinkage. Vertices within `pin_band` of the original surface (`surface`) are projected onto it at every step (parry3d BVH), so mounting faces and holes keep their CAD geometry. A vertex is pinned only when its projection direction agrees with the normal of a freely smoothed copy (|cos| > 0.6); pinning both sides of an edge would fold the triangles across it. If smoothing still creates crossings, the offending triangles' vertices are unpinned, and after a few rounds the iteration count is halved.
4. **Decimation** (`decimate`): quadric-error edge collapses (QEM) to the triangle budget, guarded by the link condition (manifoldness), valence ≥ 3, normal flips (> 60°), triangle quality and the allowed surface error. Pinned vertices do not move.
5. **Validation** (`Mesh::validate`): every guarantee is checked, not assumed; the report is `Report::mesh`.

```rust
use topopt_mesh::{Field, Options, Projector, extract};

// Samples at the optimizer's cell centres, with a ring of zeros around
let field = Field::from_cells([nx, ny, nz], h, grid_origin, &densities, 1);
let anchors = field.cell_mask([nx, ny, nz], 1, &pads);            // anchor cells
let surface = Projector::new(&body_positions, &body_triangles)?;  // the original surface (optional)
let (mesh, report) = extract(field, Some(&anchors), Some(&surface), &Options::default())?;
assert!(report.mesh.ok());
std::fs::write("part.stl", mesh.to_stl("topopt"))?;
```

## `Options`

| Field | Default | Meaning |
|---|---|---|
| `iso` | 0.5 | Surface threshold |
| `nudge` | 0.02 | Values closer than this to `iso` are moved off it; also the smallest edge fraction of a contour vertex. Smaller values put vertices closer than f32 resolution in STL (≈ 8 µm at 100 mm) and CAD welds them |
| `blur_passes` | 0 | Blur before contouring. Off: on thin convex members it moves the surface inwards (a torus loses 7 % of its volume, against 1.8 % without) |
| `remove_islands` | true | Remove components without an anchor |
| `fill_cavities` | true | Fill closed cavities |
| `smooth_iterations` | 10 | Taubin iterations |
| `lambda`, `mu` | 0.5, −0.53 | Taubin coefficients |
| `pin_band` | 1.0 | Width of the pinning band at the original surface, in sample spacings |
| `target_triangles` | `Some(40 000)` | Triangle budget (`None`: no decimation) |
| `max_decimation_error` | 0.25 | Largest RMS surface change of decimation, in sample spacings |
| `check_self_intersections` | true | Look for self-intersections after smoothing and in `validate` |

## `Report`

`cleanup` (components and samples removed or filled), `contour_triangles`, `contour_volume`, `smoothed_volume`, `pinned_vertices`, `smooth_iterations` (smoothing iterations kept), `decimated_from`, `collapses`, `mesh: MeshReport`.

## `Mesh` and `MeshReport`

```rust
pub struct Mesh { pub positions: Vec<[f64; 3]>, pub triangles: Vec<[u32; 3]> }
```

| Method | Effect |
|---|---|
| `validate(self_intersections) -> MeshReport` | Full check. `ok()` holds when there are triangles; no boundary, non-manifold or misoriented edges; no pinch vertices, degenerate triangles, self-intersections, or vertices that coincide once rounded to f32; and the volume is positive |
| `volume()`, `area()`, `mean_edge()` | Signed volume (positive with outward normals), area, mean edge length |
| `components()` | Connected components by shared edges |
| `compact()` | Drop unused vertices |
| `face_normal(t)` | A triangle's normal |
| `to_stl(header)` | Binary STL |
| `to_ply(comment)` | Binary little-endian PLY (indexed, about half the size of STL) |

`MeshReport`: `vertices`, `triangles`, `edges`, `boundary_edges`, `nonmanifold_edges`, `misoriented_edges`, `nonmanifold_vertices`, `unused_vertices`, `degenerate_triangles`, `components`, `euler`, `genus` (through-holes), `volume`, `area`, `self_intersections` (`None` when not checked), `f32_collisions` (distinct vertices equal in f32, which an STL reader would weld).

## Reading (`io`)

- `read_stl(bytes, scale)`: binary or ASCII STL. Corners that are equal bit for bit (f32) are welded into one vertex, so a closed STL gives a closed indexed mesh; triangles that collapse in welding are dropped. `scale` multiplies coordinates (1000 for a file in metres).
- `read_ply(bytes)`: the binary PLY that `to_ply` writes; other PLY layouts are refused.

## Other modules

| Module | Public items |
|---|---|
| `field` | `Field { n, origin, h, values }`: `from_cells`, `cell_mask`, `blur`, `nudge`, `zero_border`, `components`; `clean_components` |
| `mt` | `marching_tetrahedra(&Field, iso) -> Mesh` |
| `smooth` | `taubin(mesh, iters, λ, μ, constrain)`, `Projector` (closest point on a surface, parry3d BVH) |
| `decimate` | `decimate(mesh, pinned, &DecimateOptions)`: `target_triangles`, `max_error`, `max_normal_turn_deg` (60°), `min_quality` (0.15) |

## Verification

| Test | Shows |
|---|---|
| Sphere from an SDF | Closed, manifold, genus 0, volume within 3 % of 4/3 πr³ |
| Torus, whole pipeline | Genus 1, one component, volume −1.8 % |
| 60 random fields, 4³…10³, raw contour and whole pipeline | Always passes `validate`, self-intersection check included; exactly one component after clean-up |
| Pinned vertices | Lie on the original surface (error < 0.05 h) |
| Component clean-up | Islands removed, cavities filled |
| STL and PLY round trips | Welding, volume kept, ASCII and scale |

End to end as well: trimesh (`tools/reference/check_stl.py`: watertight, consistent winding, volume) and SolidCraft (`file.insert_mesh` → `closed: true` → `mesh.convert` to a solid; a 30 000-triangle cantilever converted in 5.1 s with its volume unchanged to 1e-10).

```sh
cargo +stable-x86_64-pc-windows-gnu test --release -p topopt-mesh

# Optimize a cantilever and write an STL: OUT.stl [NX NY NZ] [H_mm] [TRIANGLES]
cargo +stable-x86_64-pc-windows-gnu run --release -p topopt-mesh --example cantilever_mesh -- cantilever.stl 48 24 24 2 30000
```

## Limitations

- **Stepped domain edges.** A pinned vertex is projected onto the nearest face of the original surface, so the line where two faces meet (a box edge, a hole's rim) stays slightly stepped. Detecting sharp edges and projecting the band along them onto the edge would fix it.
- Surface detail is bounded by the field spacing: a member one or two samples thick comes out lumpy or breaks, which is why the optimizer bounds member size from below with its filter (`min_member`).
- `mesh.convert` in SolidCraft makes one planar face per triangle, so keep budgets around 20–40k triangles.
