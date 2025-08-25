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
   - `--epsilon` sets custom epsilon (LJ strength) for each bead pair.
   - Optionally add `--seed 12345` for reproducible randomization.

3. **Output structure:**  
   - Folder with all required files:
     - Randomized polymer coordinate files (`X.gro`, `P.gro`, `C.gro`)
     - Topology (`afm_system.top`)
     - All required `.itp` files (with custom LJ params)
     - Ready-to-run MDP files (`EM.mdp`, `EQ.mdp`, `production.mdp`)
     - `run.sh` script for GROMACS

4. **Run your simulation:**  
   ```bash
   cd <your_output_folder>
   bash run.sh
   ```

## Notes

- The builder script currently copies template .gro files—plug in your own randomization for full system assembly.
- The script will log the random seed and parameters in `afm_build.log`.
- Update polymer counts with `--nxylo`, `--npctn`, `--ncell`.

## Advanced

- Edit `afm_build_sweep.py` to add more control, config file support, or extend with new bead types.
- All code is pure Python 3 and requires only the standard library.

## Citation

Please cite the Sudowoodo FF if using in publications.