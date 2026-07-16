# Sudowoodo AFM Builder

This repository contains a streamlined GROMACS system builder for AFM-based cell wall simulations.

## Quick Start

1. **Prepare template files:**  
   Place the following in your repo root (or update the script to point elsewhere):  
   - `X.gro` (Xyloglucan chain)
   - `P.gro` (Pectin chain)
   - `C.gro` (Cellulose chain)
   - `toppar_custom/sudowoodo_base.itp`
   - `toppar_custom/sudowoodo_xyloglucan.itp`
   - `toppar_custom/sudowoodo_pectin.itp`
   - `toppar_custom/sudowoodo_cellulose.itp`

2. **Run the builder script:**  
   ```bash
   python afm_build_sweep.py --out run_$(date +%s) --epsilon CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4
   ```

   - `--out` specifies the output folder.
   - `--epsilon` sets custom epsilon (LJ strength) for the base C/X/P bead pairs. `PP` in `--epsilon` sets the `P/P` (neutral pectin) pair.
   - Optionally add `--seed 12345` for reproducible randomization.
   - Optionally add `--multilayer` to generate a 4-layer fiber system (see below).
   - Optionally add `--deform "0 0 0.0001 0 0 0"` to write a deformation tensor into `production.mdp` for z-axis loading.

### Per-Bead Pectin Variant Builder

For per-bead pectin epsilon assignments in `0.1` steps (`0.1` to `5.0`), use:

```bash
python build_sweep.py --out run_$(date +%s)
```

Optional flags:
- `--chains` to set the number of pectin chains to assign (default `1`)
- `--seed` for reproducible assignments
- `--pectin-itp-template` to override the pectin template path (default `toppar_custom/sudowoodo_pectin.itp`)

This writes:
- `sudowoodo_base.itp` (per-bead pectin atomtypes and nonbonded parameters)
- `sudowoodo_pectin.itp` (copied template)
- `pectin_assignment_report.txt` (sorted assignment report)

3. **Output structure:**  
   - Folder with all required files:
     - Randomized polymer coordinate files (`X.gro`, `P.gro`, `C.gro`)
     - Topology (`afm_system.top`)
     - All required `.itp` files (with custom LJ params)
     - Ready-to-run MDP files (`EM.mdp`, `EQ.mdp`, `production.mdp`)
     - `run.sh` script for GROMACS

### Monomeric Pectin Sweep

For a sweep with 100 identical unbonded pectin beads and a pectin-pectin epsilon range from `0.1` to `5.0` in `0.1` increments, use:

```bash
python build_pectin_monomer_sweep.py --out pectin_monomer_sweep
```

This creates one case directory per epsilon value, each with:
- `100` single-bead pectin molecules (no bonds, so they are not fibers)
- a `production.mdp` configured for a `100 ns` production run
- a case-specific `P-P` epsilon in `toppar_custom/sudowoodo_base.itp`

After the runs finish, you can analyze the average nearest-neighbor count for each case with:

```bash
cd pectin_monomer_sweep
python ../analyze_monomer_neighbors.py
```

The analysis script lives in the repository root and scans the current working directory for `pp_eps_*` folders.

By default it:
- uses a `5.0 nm` cutoff
- analyzes the last `25%` of frames from `production.xtc`
- averages the number of neighbors within the cutoff for each bead
- saves `nearest_neighbors_vs_epsilon.png`

Adjust the analysis constants at the top of the script to change the cutoff, frame fraction, or file names.

4. **Run your simulation:**  
   ```bash
   cd <your_output_folder>
   bash run.sh
   ```

## Deformation Runs

The builder supports adding a `deform` line to `production.mdp` with the `--deform` flag.

Example for uniaxial deformation along the z axis:

```bash
python afm_build_sweep.py --out run_$(date +%s) --epsilon CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4 --deform "0 0 0.0001 0 0 0"
```

The deform tensor follows GROMACS ordering:

```text
xx yy zz xy xz yz
```

For the production run template, deformation is intended along **z**. The generated `production.mdp` uses anisotropic Parrinello-Rahman pressure coupling with the z axis excluded from barostat scaling:

```ini
Pcoupl                   = parrinello-rahman
Pcoupltype               = anisotropic
tau_p                    = 12.0
ref_p                    = 1 1 1 0 0 0
compressibility          = 3e-4 3e-4 0 0 0 0
```

This means:
- `deform` drives the z-dimension of the box
- `x` and `y` remain pressure-coupled and can relax
- `z` is not barostatted, so the barostat does not fight the imposed strain along the loading axis

## Multi-Layer Mode

The builder supports creating a 4-layer fiber system using the `--multilayer` flag:

```bash
# For afm_build_sweep.py
python afm_build_sweep.py --out run_$(date +%s) --epsilon CC=1.0,CX=0.8,CP=0.7,XX=0.6,XP=0.5,PP=0.4 --multilayer

# Or directly with build_system.py
python build_system.py --seed 12345 --multilayer
```

When `--multilayer` is set:
- Four layers of fibers are stacked together with no spacing between layers
- Each layer contains the same types and numbers of fibers (pectin, cellulose, xyloglucan)
- Layers 1 and 3 have 0° rotation, while layers 2 and 4 are rotated 180° around the Z-axis
- The simulation box Z-dimension is automatically expanded by 4x to accommodate all layers

## Notes

- The builder script currently copies template .gro files—plug in your own randomization for full system assembly.
- The script will log the random seed and parameters in `afm_build.log`.
- Update polymer counts with `--nxylo`, `--npctn`, `--ncell`.

## Advanced

- Edit `build_sweep.py` to add more control, config file support, or extend with new bead types.
- `build_sweep.py` uses `Union[...]` type-hint syntax where needed for compatibility with Python versions prior to 3.10.
- All code is pure Python 3 and requires only the standard library (plus numpy, scipy, tqdm for build_system.py).

## Citation

Please cite the Sudowoodo FF if using in publications.
