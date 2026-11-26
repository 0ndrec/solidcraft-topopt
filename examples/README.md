# Example studies

Specs for `topopt-solver` ([format](../docs/topopt-solver.md#4-the-spec-format)). Solve one with:

```sh
topopt-solver check examples/cantilever.json      # the grid and boundary conditions, in a second
topopt-solver solve examples/cantilever.json --out examples/out
```

| Spec | Study | Result |
|---|---|---|
| [cantilever.json](cantilever.json) | A 160×60×60 mm aluminium block clamped at x = 0; 1000 N down along the free bottom edge; keep 15 % | 59 584 elements, 200 iterations, ≈ 170 s on 3 threads; 13.9 % of the volume, 216 g — [picture](../docs/images/cantilever.png) |
| [bridge.json](bridge.json) | A 240×40×60 mm steel span on two bottom supports (one free to slide along X), 5000 N spread over the deck; keep 12 % | 70 800 elements, 200 iterations, ≈ 210 s on 6 threads; 8.8 % of the volume, 398 g — [picture](../docs/images/bridge.png) |

`topopt-solver example` prints a smaller cantilever with every field at its default: a starting point for your own specs. The pictures were rendered with `scripts/render_mesh.py`.
