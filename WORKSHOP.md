# TopOpt — topology optimization in SolidCraft

Pick the body material may fill, the faces it is held by and the faces it is pushed on, and say how much of it the part may keep. TopOpt finds the stiffest shape and inserts it into your model, as a mesh or a solid.

- **A panel guides you:** Design body, Add support, Add load (several load cases), Keep surface for faces that must stay exactly as modelled. Then Optimize, and Create body.
- **Every support stays joined.** Without help, the optimizer feeds the load into the nearest supports and may leave a far bolt hole with no material. TopOpt then optimizes again, with each such support made to hold the load, so every mounting point belongs to the part.
- **Checked results:** the result is analysed once more under every load case. The Run tab shows its mass, largest deflection (mm) and von Mises stress (MPa); compare the stress with your material's strength.
- **Clean meshes:** smooth, closed and manifold, with the CAD model's own faces where you keep surfaces, and sharp edges kept. Create body turns the mesh into a solid.
- **Materials:** aluminium, steel, titanium, PLA, PETG, ABS, nylon, or your own.

## Faster with TopOpt Solver (optional)

The plugin solves on its own, on one core of the plugin sandbox. For bigger studies, download TopOpt Solver (all CPU cores) or TopOpt CUDA Solver (NVIDIA GPU) from the [GitHub release](https://github.com/0ndrec/solidcraft-topopt/releases). Start it on the folder you granted the plugin as `results`, and Optimize hands the run to it. On an entry-level GPU, a 3-inch FPV frame took 18 minutes and went from 38.2 g to 11.8 g.

## Permissions

- **files:write** — the `results` folder you choose, where the plugin and TopOpt Solver exchange jobs and results.
- **clock** — the time limit and the progress shown while it runs.
- **commands** — read the design body's faces and mesh, ask which face a point lies on, measure bodies, insert the result as a mesh and convert it to a solid, and keep the study with the design.

No network access.

## Documentation

- [User guide](https://github.com/0ndrec/solidcraft-topopt/blob/main/docs/plugin.md)
- [Choosing the settings](https://github.com/0ndrec/solidcraft-topopt/blob/main/docs/settings-guide.md), measured on one bracket
- [FPV frame example](https://github.com/0ndrec/solidcraft-topopt/blob/main/examples/fpv-frame/README.md)

Licence: MIT OR Apache-2.0.
