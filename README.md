# attenupy

Design, build, run and harvest **penEasy/PENELOPE** radiation-shielding studies
of polymer composites, from "I want 300 random 1–3-constituent mixtures per
beam" to a CSV of transmission ratios.

You supply two executables built for your machine: `penEasy.x` and
PENELOPE's `material.exe` (in its `pendbase/`, next to `pdfiles/`).
attenupy compiles nothing and generates everything else: folders,
`penEasy.in`, the `.geo` geometry, the `.spc` spectrum, `mat.mat` and
`Aire.mat`. It uses only the Python standard library (Python ≥ 3.8).

```bash
pip install -e .                                  # adds the `attenupy` command
export ATTENUPY_PENEASY=/path/to/penEasy.x
export ATTENUPY_MATERIAL=/path/to/pendbase/material.exe
attenupy doctor                                  # checks both, incl. OS/CPU mismatch
```

## Quick start

```python
import attenupy as atp
atp.log_to_console()   # print progress ("Prepared N folders...", etc.) as the study runs

# A design is just a list of simulations to create; RandomSampler draws them
# randomly (n per beam), ThicknessSweep sweeps fixed mixtures over thicknesses,
# and the two combine with `+` into one design (see "Designs" below for more).
design = atp.RandomSampler(
    n=300,                                         # simulations to draw, per beam
    beams=["RQR6_080", "RQR8_100", "x120kV"],      # or "80", "100kV", ...
    constituents=["Bi2O3", "BaSO4", "BN", "B4C", "WO3", "Sb", "Sn", "Bi", "W", "Cu", "B", "Ba"],
    n_constituents=[1, 2, 3],                      # 1-3 constituents per mixture; or {2: 0.7, 3: 0.3}
    polymer_fractions=[0.20, 0.30, 0.40],          # candidate polymer (matrix) weight fractions
    thicknesses_mm=atp.thickness_range(0.15, 1.5, 0.15),  # 0.15, 0.30, ..., 1.50 mm
    seed=0,                                        # reproducible random draw
)
held_out = atp.ThicknessSweep(atp.Mixture({"Bi2O3": 0.40, "BaSO4": 0.20}, polymer_fraction=0.40),
                             atp.thickness_range(0.15, 1.5, 0.15), beams=["RQR6_080"])

# A Study ties a design to a directory; with no toolchain= it reads penEasy.x
# and material.exe from the ATTENUPY_PENEASY / ATTENUPY_MATERIAL env vars set above.
study = atp.Study("runs/paper")
study.prepare(design + held_out)      # `+` concatenates the designs; writes every folder's inputs
study.run(max_concurrent=10)          # resumable; safe with several jobs at once
study.collect(output="results.csv")   # molecule_i / ratio_i columns by default
```

Or without Python, from a design file ([`examples/design.toml`](examples/design.toml)):

```bash
attenupy prepare runs/paper --design examples/design.toml --dry-run   # preview
attenupy prepare runs/paper --design examples/design.toml
attenupy status  runs/paper
attenupy run     runs/paper --max-concurrent 10
attenupy collect runs/paper -o results.csv            # --tidy for w_<name> columns
```

## Examples

New here? Read them in this order: [`03_gd_pb_sweep.py`](examples/03_gd_pb_sweep.py)
is the shortest complete example (design → prepare → run → collect in ~20
lines of actual code); [`01_paper_dataset.py`](examples/01_paper_dataset.py)
shows a full study with a random design plus held-out curves; and
[`02_extending.py`](examples/02_extending.py) shows how to plug in your own
spectrum, splitter, sampler and result column. Every script is heavily
commented and runnable as-is — `atp.log_to_console()` at the top of each
prints its progress so you can see what it's doing.

| file | what it does |
|---|---|
| [`01_paper_dataset.py`](examples/01_paper_dataset.py) | the paper's training design plus held-out sweeps, prepared, run and collected |
| [`02_extending.py`](examples/02_extending.py) | custom constituent, beam/spectrum, splitter, sampler and result column (no executables needed) |
| [`03_gd_pb_sweep.py`](examples/03_gd_pb_sweep.py) | Gd and Pb (70 wt%) in 30 wt% C2H4 at 100 kV, 10 thicknesses from 0.15 to 1.5 mm: 20 simulations |
| [`design.toml`](examples/design.toml), [`gd_pb_sweep.toml`](examples/gd_pb_sweep.toml) | the same kind of designs as files for `attenupy prepare` |

## Concepts

| | |
|---|---|
| `Mixture` | shielding constituents + polymer, weight fractions summing to 1 |
| `SimulationSpec` | one run: beam + mixture + thickness (mm) + penEasy seeds |
| sampler | anything that yields specs; samplers combine with `+` |
| `Study` | one directory tree; turns specs into folders, runs them, collects rows |

Every folder stores its exact spec in `spec.json`, so later stages never have
to reverse-engineer parameters from folder names. The layout is
`<root>/<beam folder>/<name>/` (e.g. `80kV/RQR6_080_2_04_oneL080_BaSO4_Bi2O3_0.20_0.40`),
plus `<root>/.attenupy/` for the material cache and provenance.

### Designs

- **`RandomSampler`**: `n` random mixtures per beam, one random thickness each.
  Mixture size, combinations (`"balanced"` spreads draws evenly over sizes and
  combinations; `"uniform"` draws independently), `require_any=["B", "Ba"]`,
  and the mass split (`splitter=`) are all configurable. Fractions are rounded
  to 2 decimals (`decimals=`) while keeping the exact total and no zero shares.
  `unique="spec"` (default) never repeats a simulation; `unique="mixture"`
  never repeats a composition (the paper's rule).
- **`ThicknessSweep`**: given mixture(s), every thickness in a list. Also
  works with random mixtures: `atp.RandomMixtures([...]).sample(6, seed=1)`.
- **`GridSampler`**: beams × mixtures × thicknesses.
- **`CSVSampler`**: specs from CSV rows. The default parser reads a wide CSV
  layout (`molecule_i`/`ratio_i` columns, see `wide_csv_row`); pass
  `parse_row=` for any other layout (it may return several specs per row).
- **`Dirichlet(alpha=1.0, min_share=0.0)`**: the flat Dirichlet of the paper;
  `Dirichlet({"B": 2.5, "Ba": 2.5}, min_share=0.1)` biases the split toward
  chosen elements while reserving a floor share for each one. `EqualSplit()`
  is also built in.

### Inputs attenupy generates

- **Geometry**: `SlabGeometry()` writes the slab/detector/air-sphere setup
  (a square slab facing a point source, an air detector behind it, all inside
  an air sphere); every dimension (`slab_half_width_cm`, `detector_z_cm`, …)
  is a parameter. Use `FileGeometry(dir)` to copy your own `.geo` files
  instead.
- **Spectra**: five beams are built in with their spectra bundled
  (`N30_030`, `RQR5_070`, `RQR6_080`, `RQR8_100`, `x120kV`).
  New beams can compute theirs with SpekPy (`pip install spekpy`):

  ```python
  rqr5 = atp.SpekpySpectrum(kvp=70, th=12, filters=[("Al", 2.8), ("Be", 3.0)])  # dk=1 keV, z=100 cm
  rqr5.hvl()                        # (2.581, 3.652) mm Al; IEC 61267 RQR5 asks for 2.58
  atp.register_beam("RQR5_spekpy", kvp=70, spectrum=rqr5, folder="70kV_spekpy")
  ```

  Each row is a SpekPy bin (lower edge in eV, unrounded fluence) and the file
  ends with `<kVp in eV>  -1`. Or use a `.spc` file
  (`spectrum="W60.spc"`), `TabulatedSpectrum(energies_eV, weights)`, or any
  object with `write(path)`. A spectrum given as point values f(E) (what most
  calculators print) should go through `TabulatedSpectrum.from_points`, which
  builds the histogram penEasy expects (manual, Appendix B). Spectrum files are
  checked with `atp.check_bigs_spectrum`, which applies penEasy's own reading
  rules (`#` comments in column 1, a negative-weight terminator, non-decreasing
  energies, at most 31999 bins) before anything is written. Only relative
  weights matter: penEasy normalises them.
- **Deck**: `PeneasyDeck()` renders `penEasy.in` for **penEasy v.2015-05-30**
  (the section versions that build expects, material list ended by `0`).
  Histories, time limit, dumps, the source cone half-angle and the
  fluence-tally range are parameters. Newer penEasy releases (e.g. 2020)
  changed section versions and end the material list with `-1`: for those,
  pass that release's `penEasy.in` to `PeneasyDeck.from_template_file()`.
- **Materials**: `mat.mat` via `material.exe` from the element weight
  fractions and the additivity-of-volumes density; `Aire.mat` from PENELOPE's
  database (#104). Materials are cached by composition, so a 10-point
  thickness sweep calls `material.exe` once, and unique materials are built
  in parallel.
- **Constituents** beyond the built-in ones:
  `atp.register_constituent("Gd2O3", density=7.41)`, or with a label that
  differs from the formula: `("PU", density=1.2, formula="C3H8N2O")`.
  Formulas may contain parentheses (`Ca10(PO4)6(OH)2`).

### I0 and the transmission ratio

`energy_ratio` = detector energy / I0, where I0 is the detector signal with
the slab replaced by air. The five bundled beams carry I0 values valid for
the default geometry, so nothing else is needed to collect `energy_ratio`
for them.

For any other beam -- one you registered yourself, or a bundled beam used
with a non-default geometry -- I0 isn't known yet, and `collect()` can't
just guess it. Compute it by adding an air-reference simulation and running
it, the same way as any other simulation:

```python
study.prepare(design)                 # your actual composite simulations
study.prepare_references()            # + one air-filled run per beam used above
study.run(max_concurrent=10)          # runs both together
study.collect(output="results.csv")   # I(t) for the composites, I0 from the air runs
```

`prepare_references(beams=None, thickness_mm=0.30)` defaults to every beam
already used in the study, at a 0.30 mm air "slab" -- the exact thickness
barely matters, since air attenuates negligibly at these energies over a few
millimetres. Pass `beams=["W60"]` to add a reference for specific beams only;
`attenupy references ROOT` does the same from the CLI.

Once a reference run for a beam has finished, `collect()` uses its measured
I0 automatically and in preference to any baked-in value -- there's nothing
extra to tell it. If neither a measured nor a baked-in I0 exists for a beam,
`collect()` prints a one-line warning and leaves that beam's `energy_ratio`
column empty rather than guessing at a number.

### Running

```bash
attenupy run runs/paper --max-concurrent 10     # run again to resume
```

Runs are tracked by small files inside each folder, so an interrupted study
picks up where it stopped, and two runners on the same study never start the
same simulation.

- A simulation is **done** when `penEasy.out` ends with penEasy's sign-off. A
  crash leaves it **failed** (retry with `--retry-failed`) instead of silently
  counting as processed.
- Claims left by dead processes on the same host are reclaimed automatically;
  `attenupy reset ROOT --state running` clears ones from other nodes.
- SIGTERM (e.g. a scheduler's time limit) stops penEasy and releases the claims.

## Extending

Every stage is duck-typed: pass any object with the right method(s) to `Study`
or use it in a design. See [`examples/02_extending.py`](examples/02_extending.py)
(runnable without executables).

| stage | interface | built in |
|---|---|---|
| design | iterable of `SimulationSpec` (subclass `Sampler` for `+`) | `RandomSampler`, `GridSampler`, `ThicknessSweep`, `CSVSampler` |
| composition | `draw(count, rng) -> [Mixture]` | `RandomMixtures` |
| mass split | `split(names, total, rng) -> [float]` | `Dirichlet`, `EqualSplit` |
| geometry | `filename(spec)`, `write(spec, path)`, optional `validate(spec)` | `SlabGeometry`, `FileGeometry` |
| spectrum | `write(path)`, optional `validate()` | `FileSpectrum`, `TabulatedSpectrum`, `SpekpySpectrum` |
| materials | `build(mixture, path, registry)`, `build_air(path)` | `PenelopeMaterials` |
| folder names | `name(spec) -> str` | `CompactNamer`, `DescriptiveNamer` |
| runner | `run(dirs) -> RunSummary` | `LocalRunner` |
| result columns | `extractor(sim, ctx) -> dict` | `WideColumns`, `SpecColumns`, `material_properties`, `EnergyRatio` |

Custom samplers become usable in design files with `@atp.register_sampler("name")`.

Material numbering is a contract between geometry and deck: 1 = slab,
2 = air detector, 3 = air enclosure.

## Safety checks

- A whole design is validated before anything is written: unknown beams or
  constituents, a slab thicker than the gap to the detector, and file names
  over penEasy's limits are all reported together.
- Two different specs mapping to the same folder name get `__2`, `__3`
  suffixes; identical specs are kept once.
- The source cone must fit on the slab (`SlabGeometry.check_beam`); otherwise
  part of the beam would reach the detector unshielded.
- penEasy seeds must be positive (0,0 and -1,-1 mean "seeds file" and
  "restart" to penEasy); an empty input file (penEasy creates one when a
  material is missing) counts as missing.
- `CompactNamer` refuses polymer fractions that aren't multiples of 0.1 (it
  would mislabel them); use `DescriptiveNamer()` for such designs.
- `collect` warns when a folder's `mat.mat` density disagrees with its spec
  by more than 1%, i.e. the simulated material is not the labelled one.
- `material.exe` runs in a private temporary directory, which avoids its
  80-character output-path limit (easy to hit with long multi-constituent
  folder names), works with read-only installs, and allows parallel builds.

## Folders without a `spec.json`

A folder made outside attenupy is recognised from its name (the
`CompactNamer` scheme), so `attenupy status`, `materials` and `collect` work
on it directly. A folder that carries its own `penEasy.x` runs with it when
no `--peneasy` is given.

## Tests

```bash
python -m unittest discover -s tests
```

The tests use a stand-in `penEasy.x`. The `material.exe` test runs when
`$ATTENUPY_MATERIAL` points to a working one; the SpekPy tests run when
`spekpy` is installed.
