#!/usr/bin/env python3
"""Sweep over the number of PR and PC beads per fiber in a full cell-wall system.

Epsilons for PR, PN, and PC are fixed by the caller.  The script iterates
over PR bead counts (0, 5, 10, …, up to --pr-max) and PC bead counts
(0, 5, 10, …, up to --pc-max) in steps of --bead-step.  Any remaining
beads in a fiber are assigned to PN so that every fiber always totals
BEADS_PER_FIBER (30) beads.  Cases where pr_count + pc_count >
BEADS_PER_FIBER are skipped.

For each (pr_count, pc_count) combination the script calls
build_system.build() which constructs the full AFM cell-wall system
(cellulose + xyloglucan + pectin) inside its own output subdirectory.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

# Allow importing siblings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import afm_build_sweep
import build_sweep
import build_system


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep over PR and PC bead counts per fiber, building the full "
            "AFM cell-wall system for each combination. "
            "PN beads fill the remainder up to BEADS_PER_FIBER (30). "
            "Epsilons for PR, PN, and PC are fixed."
        )
    )
    parser.add_argument("--out", type=Path, required=True, help="Output folder")
    parser.add_argument(
        "--pr-epsilon", type=float, default=build_sweep.DEFAULT_EPSILON_BY_TYPE["PR"],
        help=f"PR bead epsilon in kJ/mol (default {build_sweep.DEFAULT_EPSILON_BY_TYPE['PR']})",
    )
    parser.add_argument(
        "--pn-epsilon", type=float, default=build_sweep.DEFAULT_EPSILON_BY_TYPE["PN"],
        help=f"PN bead epsilon in kJ/mol (default {build_sweep.DEFAULT_EPSILON_BY_TYPE['PN']})",
    )
    parser.add_argument(
        "--pc-epsilon", type=float, default=build_sweep.DEFAULT_EPSILON_BY_TYPE["PC"],
        help=f"PC bead epsilon in kJ/mol (default {build_sweep.DEFAULT_EPSILON_BY_TYPE['PC']})",
    )
    parser.add_argument(
        "--pr-max", type=int, default=25,
        help="Maximum number of PR beads per fiber (default 25)",
    )
    parser.add_argument(
        "--pc-max", type=int, default=25,
        help="Maximum number of PC beads per fiber (default 25)",
    )
    parser.add_argument(
        "--bead-step", type=int, default=5,
        help="Increment step for PR and PC bead counts (default 5)",
    )
    parser.add_argument("--seed", type=int, default=None, help="Random seed (default: current time)")
    parser.add_argument(
        "--multilayer", action="store_true",
        help="Generate 4-layer fiber system (passed through to build_system.build)",
    )
    parser.add_argument("--gmx", type=str, default="gmx", help="GROMACS launcher (default: gmx)")
    parser.add_argument("--ntomp", type=int, default=24, help="Number of OpenMP threads (default: 24)")
    parser.add_argument("--ntmpi", type=int, default=1, help="Number of MPI ranks (default: 1)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Sweep helpers
# ---------------------------------------------------------------------------

def iter_bead_counts(pr_max: int, pc_max: int, step: int):
    """Yield (pr_count, pc_count) pairs that fit within BEADS_PER_FIBER."""
    if step <= 0:
        raise ValueError("bead-step must be positive")
    for pr in range(0, pr_max + 1, step):
        for pc in range(0, pc_max + 1, step):
            if pr + pc <= build_sweep.BEADS_PER_FIBER:
                yield pr, pc


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    base_seed = args.seed if args.seed is not None else int(time.time())

    cases = list(iter_bead_counts(args.pr_max, args.pc_max, args.bead_step))
    print(f"[INFO] Sweeping {len(cases)} (PR, PC) combinations")
    print(f"[INFO] Epsilons — PR: {args.pr_epsilon}  PN: {args.pn_epsilon}  PC: {args.pc_epsilon}")

    orig_dir = os.getcwd()

    for pr_count, pc_count in cases:
        pn_count = build_sweep.BEADS_PER_FIBER - pr_count - pc_count
        case_name = f"pr{pr_count:02d}_pc{pc_count:02d}_pn{pn_count:02d}"
        case_dir = args.out / case_name
        case_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"[CASE] {case_name}  (PR={pr_count}, PC={pc_count}, PN={pn_count})")
        print(f"{'='*60}")

        # build_system.build() writes files into cwd, so we switch to the case dir
        os.chdir(case_dir)
        try:
            build_system.build(
                seed=base_seed,
                multilayer=args.multilayer,
                pr_epsilon=args.pr_epsilon,
                pn_epsilon=args.pn_epsilon,
                pc_epsilon=args.pc_epsilon,
                pc_per_fiber=pc_count,
                pr_per_fiber=pr_count,
                src_dir=orig_dir,
            )
        finally:
            os.chdir(orig_dir)

        afm_build_sweep.write_mdp_files(args, case_dir)
        afm_build_sweep.write_run_sh(args, case_dir)
        print(f"[ok] finished {case_dir}")

    print(f"\n[DONE] All {len(cases)} cases written to {args.out}")


if __name__ == "__main__":
    main()
