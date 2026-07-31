#!/usr/bin/env python3
"""
Sweep over a grid of PR-epsilon × PC-epsilon values.

For each (PR epsilon, PC epsilon) pair this script:
  1. Calls afm_build_sweep.py to build the run directory.
  2. Writes a SLURM submit script (submit.sh) into that directory.
  3. Submits the job with sbatch (unless --dry-run is given).

Directory names follow the pattern:
  <out-root>/pr<PR_EPS>_pc<PC_EPS>
  e.g.  sweep/pr0.4_pc4.8

Usage examples
--------------
# Sweep PR 0.1–1.0 (step 0.1) × PC 2.0–6.0 (step 1.0), dry run:
  python run_pectin_epsilon_sweep.py \\
      --pr-start 0.1 --pr-stop 1.0 --pr-step 0.1 \\
      --pc-start 2.0 --pc-stop 6.0 --pc-step 1.0 \\
      --out-root ./sweep \\
      --dry-run

# Real run with a fixed seed and custom core epsilons:
  python run_pectin_epsilon_sweep.py \\
      --pr-start 0.2 --pr-stop 0.8 --pr-step 0.2 \\
      --pc-start 3.0 --pc-stop 5.0 --pc-step 1.0 \\
      --out-root ./sweep \\
      --seed 42 \\
      --epsilon CC=2.5,CX=2.5,CP=2.5,XX=2.5,XP=2.5,PP=2.5
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import textwrap
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import build_sweep

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_EPSILON = "CC=2.5,CX=2.5,CP=2.5,XX=2.5,XP=2.5,PP=2.5"
DEFAULT_PN_EPSILON = build_sweep.DEFAULT_EPSILON_BY_TYPE["PN"]

SLURM_TEMPLATE = textwrap.dedent("""\
    #!/bin/bash
    #SBATCH --job-name={job_name}
    #SBATCH --time=72:00:00
    #SBATCH --partition=astro2_long
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=24
    #SBATCH --threads-per-core=1
    #SBATCH --chdir={run_dir}

    module purge
    module load astro

    export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
    #export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

    source /software/astro/gromacs/2023.3-cpu-no-mpi/bin/GMXRC

    bash run.sh
""")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _iter_range(start: str, stop: str, step: str) -> list[Decimal]:
    """Return inclusive list of Decimal values from start to stop by step."""
    s, e, d = Decimal(start), Decimal(stop), Decimal(step)
    if d <= 0:
        raise ValueError("step must be positive")
    if e < s:
        raise ValueError("stop must be >= start")
    values: list[Decimal] = []
    v = s
    while v <= e + d / 2:  # half-step tolerance for float endpoints
        values.append(v)
        v += d
    return values


def _fmt(val: Decimal) -> str:
    """Format a Decimal for directory names — drop trailing zeros."""
    return str(val.normalize())


def _dir_name(pr_eps: Decimal, pc_eps: Decimal) -> str:
    return f"pr{_fmt(pr_eps)}_pc{_fmt(pc_eps)}"


def _job_name(pr_eps: Decimal, pc_eps: Decimal) -> str:
    return f"sw_pr{_fmt(pr_eps)}_pc{_fmt(pc_eps)}"


def _write_slurm_script(run_dir: Path, pr_eps: Decimal, pc_eps: Decimal, dry_run: bool = False) -> Path:
    script = SLURM_TEMPLATE.format(
        job_name=_job_name(pr_eps, pc_eps),
        run_dir=str(run_dir.resolve()),
    )
    submit_path = run_dir / "submit.sh"
    if dry_run:
        print(f"    (dry-run) would write: {submit_path}")
    else:
        submit_path.write_text(script)
        submit_path.chmod(0o755)
    return submit_path


def _build_variant(
    run_dir: Path,
    pr_eps: Decimal,
    pc_eps: Decimal,
    pn_eps: float,
    epsilon: str,
    seed: int | None,
    extra_args: list[str],
    dry_run: bool,
) -> bool:
    """Call afm_build_sweep.py to populate run_dir. Returns True on success."""
    builder = Path(__file__).parent / "afm_build_sweep.py"
    cmd = [
        sys.executable, str(builder),
        "--out", str(run_dir),
        "--epsilon", epsilon,
        "--pr-epsilon", str(pr_eps),
        "--pn-epsilon", str(pn_eps),
        "--pc-epsilon", str(pc_eps),
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    cmd.extend(extra_args)

    print(f"  [build] {run_dir.name}")
    if dry_run:
        print(f"    (dry-run) would run: {' '.join(cmd)}")
        return True

    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"  [ERROR] build failed for {run_dir.name} (exit {result.returncode})")
        return False
    return True


def _submit_job(submit_path: Path, dry_run: bool) -> None:
    print(f"  [submit] {submit_path.parent.name}")
    if dry_run:
        print(f"    (dry-run) would run: sbatch {submit_path}")
        return
    result = subprocess.run(["sbatch", str(submit_path)])
    if result.returncode != 0:
        print(f"  [ERROR] sbatch failed for {submit_path.parent.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sweep PR-epsilon × PC-epsilon, build each variant, submit SLURM jobs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # PR range
    p.add_argument("--pr-start", required=True, help="PR epsilon start value")
    p.add_argument("--pr-stop",  required=True, help="PR epsilon stop value (inclusive)")
    p.add_argument("--pr-step",  required=True, help="PR epsilon step size")
    # PC range
    p.add_argument("--pc-start", required=True, help="PC epsilon start value")
    p.add_argument("--pc-stop",  required=True, help="PC epsilon stop value (inclusive)")
    p.add_argument("--pc-step",  required=True, help="PC epsilon step size")
    # Fixed PN epsilon
    p.add_argument("--pn-epsilon", type=float, default=DEFAULT_PN_EPSILON,
                   help="Fixed PN (neutral) bead epsilon")
    # Output
    p.add_argument("--out-root", type=Path, default=Path("sweep"),
                   help="Root directory; one subdirectory is created per case")
    # Build options forwarded to afm_build_sweep.py
    p.add_argument("--epsilon", default=DEFAULT_EPSILON,
                   help="Core LJ epsilon map passed to afm_build_sweep.py")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed forwarded to afm_build_sweep.py (same for all cases)")
    p.add_argument("--multilayer", action="store_true",
                   help="Pass --multilayer to afm_build_sweep.py")
    p.add_argument("--ktheta", default=None,
                   help="Pass --ktheta value to afm_build_sweep.py")
    # Control
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without executing them")
    p.add_argument("--build-only", action="store_true",
                   help="Build directories but do not submit SLURM jobs")
    p.add_argument("--submit-only", action="store_true",
                   help="Submit existing directories without rebuilding")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    pr_values = _iter_range(args.pr_start, args.pr_stop, args.pr_step)
    pc_values = _iter_range(args.pc_start, args.pc_stop, args.pc_step)
    total = len(pr_values) * len(pc_values)
    print(f"Sweep: {len(pr_values)} PR values × {len(pc_values)} PC values = {total} cases")
    print(f"Output root: {args.out_root.resolve()}")
    if args.dry_run:
        print("*** DRY RUN — no files will be written or jobs submitted ***")

    args.out_root.mkdir(parents=True, exist_ok=True)

    # Extra flags forwarded verbatim to afm_build_sweep.py
    extra: list[str] = []
    if args.multilayer:
        extra.append("--multilayer")
    if args.ktheta:
        extra.extend(["--ktheta", args.ktheta])

    run_dirs: list[Path] = []

    # --- Build phase ---
    if not args.submit_only:
        print("\n=== Building run directories ===")
        failed: list[str] = []
        for pr_eps in pr_values:
            for pc_eps in pc_values:
                run_dir = args.out_root / _dir_name(pr_eps, pc_eps)
                ok = _build_variant(
                    run_dir, pr_eps, pc_eps, args.pn_epsilon,
                    args.epsilon, args.seed, extra, args.dry_run,
                )
                if ok:
                    run_dirs.append(run_dir)
                else:
                    failed.append(run_dir.name)
        if failed:
            print(f"\n[WARNING] {len(failed)} build(s) failed: {', '.join(failed)}")
    else:
        # Collect existing directories matching the naming pattern
        for pr_eps in pr_values:
            for pc_eps in pc_values:
                run_dir = args.out_root / _dir_name(pr_eps, pc_eps)
                if run_dir.is_dir():
                    run_dirs.append(run_dir)
                else:
                    print(f"  [skip] {run_dir.name} not found, skipping submission")

    if args.build_only:
        print(f"\nBuild-only mode: {len(run_dirs)} director(ies) ready.")
        return

    # --- Submit phase ---
    print(f"\n=== Writing submit scripts and submitting {len(run_dirs)} job(s) ===")
    for run_dir in run_dirs:
        # Recover epsilon values from directory name for labelling
        name = run_dir.name
        parts = name.split("_")
        pr_eps = Decimal(parts[0][2:])  # strip leading "pr"
        pc_eps = Decimal(parts[1][2:])  # strip leading "pc"
        submit_path = _write_slurm_script(run_dir, pr_eps, pc_eps, args.dry_run)
        _submit_job(submit_path, args.dry_run)

    print("\nDone.")


if __name__ == "__main__":
    main()
