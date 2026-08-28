# Sudowoodo AFM Builder

A GROMACS system builder for AFM-based plant cell-wall simulations. Assembles coarse-grained cellulose, xyloglucan, and pectin chains with configurable Lennard-Jones parameters, then generates all files needed to run energy minimisation, equilibration, and production MD.

---

## Repository overview

| Script | Purpose |
|---|---|
| `build_system.py` | Core system builder — places chains, writes GRO/TOP/ITP |
| `afm_build_sweep.py` | Single-case builder with custom C/X/P epsilon mapping and ktheta control |
| `build_sweep.py` | Pectin variant builder — per-fiber PR/PN/PC bead assignment |
| `build_pectin_monomer_sweep.py` | Monomeric pectin epsilon sweep (unbonded single-bead system) |
| `build_bead_count_sweep.py` | **Bead-count sweep** — full system per (PR, PC) count combination |
| `run_pectin_epsilon_sweep.py` | PR-epsilon × PC-epsilon grid sweep with SLURM submission |
| `sweep_eps.py` | Uniform epsilon-scale sweep from an existing run directory |
| `analyze_monomer_neighbors.py` | Post-processing: nearest-neighbour analysis across monomer-sweep cases |

---

## Prerequisites

```
Python ≥ 3.8
numpy, scipy, tqdm
GROMACS (gmx in PATH, or pass --gmx /path/to/gmx)
```

Install Python dependencies:

```bash
pip install numpy scipy tqdm
```

---

## Template files required

Place the following files in the directory from which you run any build script:

```
C.gro                               # Cellulose chain template
X.gro                               # Xyloglucan chain template
P.gro                               # Pectin chain template
toppar_custom/sudowoodo_base.itp
toppar_custom/sudowoodo_cellulose.itp
toppar_custom/sudowoodo_xyloglucan.itp
toppar_custom/sudowoodo_pectin.itp
```

---

## Workflows

### 1. Single system with custom LJ parameters — `afm_build_sweep.py`

Builds one run directory. You set the epsilon for each coarse-grained bead-pair (C/X/P) and optionally the angular force constants.

```bash
python afm_build_sweep.py \
  --out run_$(date +%s) \
  --epsilon CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4
```

Optional flags:

| Flag | Description |
|---|---|
| `--seed INT` | Random seed for chain placement |
| `--ktheta "p,c,x"` | Angular force constants for pectin, cellulose, xyloglucan (use empty field to keep default, e.g. `",150,"`) |
| `--multilayer` | Stack 4 layers with alternating 180° rotation |
| `--deform "0 0 e 0 0 0"` | Add a deformation tensor to `production.mdp` for z-axis loading |
| `--pr-epsilon FLOAT` | Repulsive pectin bead self-epsilon |
| `--pn-epsilon FLOAT` | Neutral pectin bead self-epsilon |
| `--pc-epsilon FLOAT` | Crosslink pectin bead self-epsilon |

#### Deformation runs

For uniaxial z-axis loading, add `--deform "0 0 0.0001 0 0 0"`. The generated `production.mdp` uses anisotropic Parrinello-Rahman coupling so that x and y relax while z is driven by the imposed strain:

```ini
Pcoupl        = parrinello-rahman
Pcoupltype    = anisotropic
ref_p         = 1 1 1 0 0 0
compressibility = 3e-4 3e-4 0 0 0 0
```

---

### 2. Direct system build — `build_system.py`

The core builder called by all other scripts. Can also be run directly:

```bash
python build_system.py \
  --seed 42 \
  --pr-epsilon 0.4 --pn-epsilon 2.2 --pc-epsilon 4.8 \
  --pr 5 --pc 2 \
  --multilayer
```

| Flag | Description |
|---|---|
| `--seed INT` | Random seed |
| `--pr-epsilon FLOAT` | PR bead self-epsilon (kJ/mol) |
| `--pn-epsilon FLOAT` | PN bead self-epsilon (kJ/mol) |
| `--pc-epsilon FLOAT` | PC bead self-epsilon (kJ/mol) |
| `--pr INT` | PR beads per fiber (default 2) |
| `--pc INT` | PC beads per fiber (default 2); remainder → PN |
| `--multilayer` | 4-layer system with Z×4 box |

Outputs are written to the current working directory:
- `afm_system.gro` / `afm_system.top`
- `toppar_custom/sudowoodo_base.itp` (pectin atomtypes + nonbonded params)
- `toppar_custom/sudowoodo_pectin_N.itp` (one per fiber)

---

### 3. Pectin bead-type variant builder — `build_sweep.py`

Generates per-fiber ITP files where each fiber has a randomised sequence of PR (repulsive), PN (neutral), and PC (crosslink) beads. Three shared atomtypes (`PctRep`, `PctNeu`, `PctXlk`) are defined with configurable epsilons.

```bash
python build_sweep.py \
  --out toppar_out \
  --chains 10 \
  --pr 5 --pc 2 \
  --pr-epsilon 0.4 --pn-epsilon 2.2 --pc-epsilon 4.8 \
  --seed 42
```

Outputs in `--out`:
- `sudowoodo_base.itp`
- `sudowoodo_pectin_1.itp … sudowoodo_pectin_N.itp`
- `pectin_assignment_report.txt`

---

### 4. Bead-count sweep — `build_bead_count_sweep.py`

**Sweeps over the number of PR and PC beads per fiber** while keeping epsilons fixed. Builds a complete full-system simulation directory for every (PR count, PC count) combination. Remaining beads are always assigned to PN so that each fiber totals 30 beads.

```bash
python build_bead_count_sweep.py \
  --out sweep_bead_counts \
  --pr-epsilon 0.4 --pn-epsilon 2.2 --pc-epsilon 4.8 \
  --pr-max 25 --pc-max 0 \
  --bead-step 5 \
  --seed 42
```

This produces one subdirectory per combination, for example with `--pr-max 25 --pc-max 0`:

```
sweep_bead_counts/pr00_pc00_pn30/   PR=0,  PC=0, PN=30
sweep_bead_counts/pr05_pc00_pn25/   PR=5,  PC=0, PN=25
sweep_bead_counts/pr10_pc00_pn20/   PR=10, PC=0, PN=20
...
sweep_bead_counts/pr25_pc00_pn05/   PR=25, PC=0, PN=5
```

Each subdirectory contains a complete GROMACS run (GRO, TOP, ITP files, MDP files, `run.sh`).

| Flag | Description |
|---|---|
| `--pr-max INT` | Maximum PR beads to sweep up to (default 25) |
| `--pc-max INT` | Maximum PC beads to sweep up to (default 25) |
| `--bead-step INT` | Increment for PR and PC counts (default 5) |
| `--pr-epsilon FLOAT` | Fixed PR bead epsilon |
| `--pn-epsilon FLOAT` | Fixed PN bead epsilon |
| `--pc-epsilon FLOAT` | Fixed PC bead epsilon |
| `--seed INT` | Random seed |
| `--multilayer` | Pass multilayer mode through to `build_system` |

Combinations where PR + PC > 30 are automatically skipped.

---

### 5. PR-epsilon × PC-epsilon grid sweep — `run_pectin_epsilon_sweep.py`

Sweeps over a grid of PR-epsilon and PC-epsilon values, building a run directory for each pair. Optionally submits each case to SLURM.

```bash
python run_pectin_epsilon_sweep.py \
  --pr-start 0.1 --pr-stop 1.0 --pr-step 0.1 \
  --pc-start 2.0 --pc-stop 6.0 --pc-step 1.0 \
  --out-root ./sweep \
  --dry-run
```

Remove `--dry-run` to submit jobs with `sbatch`.

---

### 6. Monomeric pectin epsilon sweep — `build_pectin_monomer_sweep.py`

Creates one case per P-P epsilon value (default 0.1–5.0, step 0.1). Each case contains 100 identical unbonded single-bead pectin molecules and a 100 ns production run. Useful for calibrating the pectin self-interaction epsilon before running fiber systems.

```bash
python build_pectin_monomer_sweep.py \
  --out pectin_monomer_sweep \
  --epsilon-start 0.1 --epsilon-stop 5.0 --epsilon-step 0.1 \
  --count 100 --prod-ns 100
```

Case directories are named `pp_eps_0.1`, `pp_eps_0.2`, etc.

---

### 7. Uniform epsilon-scale sweep — `sweep_eps.py`

Copies an existing run directory and scales all epsilons in `sudowoodo_base.itp` by a set of factors. Useful for sensitivity analysis around a known good configuration.

```bash
python sweep_eps.py \
  --scales 0.5 1.0 1.5 2.0 \
  --template-top afm_system.top \
  --template-gro afm_system.gro \
  --toppar-dir toppar_custom \
  --out sweep_eps
```

---

## Post-processing

### Nearest-neighbour analysis — `analyze_monomer_neighbors.py`

Analyses the monomeric pectin sweep output. For each `pp_eps_*` case directory it reads the last 25 % of frames from `production.xtc` and computes the average number of neighbours within a 1.5 nm cutoff per bead.

```bash
cd pectin_monomer_sweep
python ../analyze_monomer_neighbors.py
```

Output: `nearest_neighbors_vs_epsilon.png`

Tuneable constants at the top of the script:

| Constant | Default | Description |
|---|---|---|
| `CUTOFF_NM` | 1.5 | Neighbour cutoff (nm) |
| `LAST_FRACTION` | 0.25 | Fraction of trajectory to analyse |
| `CASE_GLOB` | `pp_eps_*` | Glob pattern for case directories |

---

## Running a simulation

Every case directory produced by the build scripts contains a `run.sh`:

```bash
cd <case_dir>
bash run.sh
```

`run.sh` sequentially runs energy minimisation (`EM`), equilibration (`EQ`), and production MD using `gmx grompp` + `gmx mdrun`.

GROMACS parallelism defaults:
- `--ntmpi 1` / `--ntomp 24` (single MPI rank, 24 OpenMP threads)

Override by passing `--ntmpi` and `--ntomp` to the build script.

---

## Multi-layer mode

Any build script that accepts `--multilayer` stacks four layers of fibers. Each layer occupies its own Z-slice of the box; layers 2 and 4 are rotated 180° around Z relative to layers 1 and 3. The simulation box Z-dimension is automatically set to 4× the single-layer Z extent.

---

## Testing

```bash
python -m unittest -v
```

For focused pectin-variant tests only:

```bash
python -m unittest -v test_pectin_variant.py
```

---

## Citation

Please cite the Sudowoodo force field if you use this repository in a publication.

