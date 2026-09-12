# PRX Fig. 6 — simulation & analysis pipeline

Vertex-model tissue simulation under tension remodeling (`MainFile.py`), plus two analysis/figure scripts
(`FigureDelayedEvents.py`, `LastSnapshot.py`), to reproduce panels as in "Pérez-Verdugo, F., & Banerjee, S. (2023). Tension remodeling regulates topological transitions in epithelial tissues. PRX life, 1(2), 023006.", Fig. 6 (b) and (c).

Column meanings for input/output .txt files are documented below by
position (0-indexed).

## Requirements

- Python 3
- `numpy`
- `matplotlib`

## Run order

```bash
python3 MainFile.py            # 1. runs the simulation, writes outputsVertices/, outputsJunctions/, outputsT1/
python3 FigureDelayedEvents.py # 2. reads outputsT1/, writes outputsT1/lifetime_unresolved.txt + FigureDelayedEvents.png
python3 LastSnapshot.py        # 3. reads outputsVertices/, outputsJunctions/ and lifetime_unresolved.txt, writes LastSnapshot.png
```

`MainFile.py` deletes and recreates `outputsVertices/`, `outputsJunctions/`,
and `outputsT1/` on every run.

`LastSnapshot.py` renders the snapshot at step index `ii=119`. The PRX figure itself used index `119`, i.e. the final snapshot. 

## Inputs — `config/`

Loaded once at the start of `MainFile.py` to build the initial tissue. 

### `config/celda.txt` 

| col | name | meaning |
|---|---|---|
| 0 | `celda` | cell index |
| 1 | `v_celda` | index of a vertex belonging to that cell |

### `config/vertices.txt` 

| col | name | meaning |
|---|---|---|
| 0 | `x` | initial x position |
| 1 | `y` | initial y position |

### `config/celulas.txt` 

| col | name | meaning |
|---|---|---|
| 0 | `area0_ec_mov` | target/reference cell area |

### `config/topology.txt` 

| col | name | meaning |
|---|---|---|
| 0 | `vtipo1` | vertex index |
| 1–6 | `adys_tipo1` | 3 neighboring-vertex pairs: `(ady1_a, ady1_b, ady2_a, ady2_b, ady3_a, ady3_b)` |
| 7–9 | `cells_tipo1` | the 3 cell indices meeting at this vertex |

## Outputs

### `outputsVertices/{d}_vertices.txt` 

Vertex positions at that timestep.

| col | name | meaning |
|---|---|---|
| 0 | `vertex_id` | vertex index |
| 1 | `x` | x position |
| 2 | `y` | y position |

### `outputsJunctions/{d}.txt` 

Cell-cell junction (edge) list at that timestep — which vertex pairs are
connected.

| col | name | meaning |
|---|---|---|
| 0 | `vertex_i` | first vertex index of the junction |
| 1 | `vertex_j` | second vertex index of the junction |

### `outputsT1/4-fold.txt`

One row per vertex merging into a transient 4-fold vertex (start of a
potential T1 transition).

| col | name | meaning |
|---|---|---|
| 0 | `t` | simulation time |
| 1 | `vi` | vertex i |
| 2 | `vj` | vertex j |
| 3 | `lij` | junction length at merge time |
| 4 | `tension` | edge tension `Tensiones[i][j]` at merge time |
| 5 | `tact` | active tension `Tact[i][j]` at merge time |
| 6 | `tension_fluct` | fluctuating tension `TensionesFluct[i][j]` at merge time |

### `outputsT1/T1_original.txt`, `outputsT1/T1_perpendicular.txt` 

One row each time a 4-fold vertex resolves back along its original axis
(`T1_original.txt`) or along the perpendicular axis (`T1_perpendicular.txt`),
i.e. a completed T1 transition of that type.

| col | name | meaning |
|---|---|---|
| 0 | `t` | resolution time |
| 1 | `vi` | vertex i (post-resolution) |
| 2 | `vj` | vertex j (post-resolution) |

### `outputsT1/lifetime_unresolved.txt` (written by `FigureDelayedEvents.py`)

For each 4-fold event that never resolves before the simulation ends
("delayed"/unresolved event).

| col | name | meaning |
|---|---|---|
| 0 | `idx` | running index |
| 1 | `lifetime` | time (minutes) the vertex stayed unresolved |
| 2 | `vertex_id` | vertex index |

### Other files written by the analysis scripts

- `FigureDelayedEvents.png` (`FigureDelayedEvents.py`) — timeline plot of
  T1/4-fold event lifetimes.
- `LastSnapshot.png` (`LastSnapshot.py`) — tissue snapshot at the chosen
  timestep, with unresolved 4-fold vertices highlighted in red.

## References

> Pérez-Verdugo, F., & Banerjee, S. (2023). Tension remodeling regulates
> topological transitions in epithelial tissues. *PRX Life*, 1(2), 023006.

## License

MIT — see [LICENSE](LICENSE).
