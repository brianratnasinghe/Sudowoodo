#!/usr/bin/env python3
"""Sweep over the number of PR and PC beads per fiber.

Epsilons for PR, PN, and PC are fixed by the caller.  The script iterates
over PR bead counts (0, 5, 10, …, up to --pr-max) and PC bead counts
(0, 5, 10, …, up to --pc-max) in steps of 5, with any remainder beads
assigned to PN so that every fiber always totals BEADS_PER_FIBER (30) beads.
Cases where pr_count + pc_count > BEADS_PER_FIBER are skipped.
"""
from __future__ import annotations

import argparse
import math
import os
import random
import shutil
import textwrap
from pathlib import Path

from build_sweep import (
    BEADS_PER_FIBER,
    DEFAULT_EPSILON_BY_TYPE,
    assign_all_chain_bead_epsilons,
    pectin_epsilon_by_type,
    write_per_bead_base_itp,
    write_per_fiber_pectin_itps,
    write_assignment_report,
)

DEFAULT_PR_MAX = 25
DEFAULT_PC_MAX = 25
DEFAULT_BEAD_STEP = 5
DEFAULT_CHAINS = 1
DEFAULT_SEED = 42
DEFAULT_BOX = (100.0, 100.0, 100.0)
DEFAULT_SPACING = 3.5
DEFAULT_PRODUCTION_NS = 100
DEFAULT_DT_PS = 0.1
PECTIN_SIGMA = 1.0
PECTIN_MASS = 26.6


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep over PR and PC bead counts per fiber. "
            "PN beads fill the remainder up to BEADS_PER_FIBER (30). "
            "Epsilons for PR, PN, and PC are fixed."
        )
    )
    parser.add_argument("--out", type=Path, required=True, help="Output folder")
    parser.add_argument(
        "--pr-epsilon", type=float, default=DEFAULT_EPSILON_BY_TYPE["PR"],
        help=f"PR bead epsilon in kJ/mol (default {DEFAULT_EPSILON_BY_TYPE['PR']})",
    )
    parser.add_argument(
        "--pn-epsilon", type=float, default=DEFAULT_EPSILON_BY_TYPE["PN"],
        help=f"PN bead epsilon in kJ/mol (default {DEFAULT_EPSILON_BY_TYPE['PN']})",
    )
    parser.add_argument(
        "--pc-epsilon", type=float, default=DEFAULT_EPSILON_BY_TYPE["PC"],
        help=f"PC bead epsilon in kJ/mol (default {DEFAULT_EPSILON_BY_TYPE['PC']})",
    )
    parser.add_argument(
        "--pr-max", type=int, default=DEFAULT_PR_MAX,
        help=f"Maximum number of PR beads per fiber (default {DEFAULT_PR_MAX})",
    )
    parser.add_argument(
        "--pc-max", type=int, default=DEFAULT_PC_MAX,
        help=f"Maximum number of PC beads per fiber (default {DEFAULT_PC_MAX})",
    )
    parser.add_argument(
        "--bead-step", type=int, default=DEFAULT_BEAD_STEP,
        help=f"Increment step for PR and PC bead counts (default {DEFAULT_BEAD_STEP})",
    )
    parser.add_argument(
        "--chains", type=int, default=DEFAULT_CHAINS,
        help=f"Number of pectin chains (default {DEFAULT_CHAINS})",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument(
        "--box", type=float, nargs=3, default=list(DEFAULT_BOX),
        metavar=("X", "Y", "Z"), help="Simulation box size in nm",
    )
    parser.add_argument(
        "--spacing", type=float, default=DEFAULT_SPACING,
        help="Grid spacing between chain COM positions in nm",
    )
    parser.add_argument(
        "--prod-ns", type=float, default=DEFAULT_PRODUCTION_NS,
        help="Production run length in ns",
    )
    parser.add_argument(
        "--dt-ps", type=float, default=DEFAULT_DT_PS,
        help="Production timestep in ps",
    )
    parser.add_argument("--gmx", type=str, default="gmx")
    parser.add_argument("--ntomp", type=int, default=24)
    parser.add_argument("--ntmpi", type=int, default=1)
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def iter_bead_counts(pr_max: int, pc_max: int, step: int):
    """Yield (pr_count, pc_count) pairs that fit within BEADS_PER_FIBER."""
    if step <= 0:
        raise ValueError("bead-step must be positive")
    pr_values = list(range(0, pr_max + 1, step))
    pc_values = list(range(0, pc_max + 1, step))
    for pr in pr_values:
        for pc in pc_values:
            if pr + pc <= BEADS_PER_FIBER:
                yield pr, pc


def production_steps(prod_ns: float, dt_ps: float) -> int:
    return int(prod_ns * 1000.0 / dt_ps)


def generate_chain_positions(chain_count: int, box: tuple, spacing: float) -> list:
    """Place chain COMs on a 2-D grid in the XY plane."""
    if chain_count <= 0:
        raise ValueError("chains must be positive")
    nx = math.ceil(math.sqrt(chain_count))
    ny = math.ceil(chain_count / nx)
    x_extent = (nx - 1) * spacing
    y_extent = (ny - 1) * spacing
    if x_extent > box[0] or y_extent > box[1]:
        raise ValueError(
            f"box {box} is too small for {chain_count} chains at spacing {spacing}"
        )
    x0 = (box[0] - x_extent) / 2.0
    y0 = (box[1] - y_extent) / 2.0
    z = box[2] / 2.0
    positions = []
    for idx in range(chain_count):
        ix = idx % nx
        iy = idx // nx
        positions.append((x0 + ix * spacing, y0 + iy * spacing, z))
    return positions


def write_gro(path: Path, assignments, chain_positions: list, box: tuple) -> None:
    """Write a GRO file with all beads from all chains."""
    from build_sweep import iter_assignments
    bead_list = list(iter_assignments(assignments))
    n_atoms = len(bead_list)
    lines = ["Pectin bead-count sweep system", f"{n_atoms:5d}"]
    chain_com = {chain_idx: pos for chain_idx, pos in enumerate(chain_positions, 1)}
    # Place each bead at the chain COM (single bead per chain for monomers,
    # or all chain beads share the same COM for simplicity in starting config).
    for bead_idx, assignment in enumerate(bead_list, 1):
        chain_idx = int(assignment["chain_index"])
        cx, cy, cz = chain_com[chain_idx]
        bead_offset = (int(assignment["bead_index"]) - 1) * 0.31  # ~bond length
        mol_name = f"Pctn_{chain_idx}"
        atom_name = f"P{assignment['bead_index']}"
        lines.append(
            f"{chain_idx:5d}{mol_name:<5}{atom_name:>5}{bead_idx:5d}"
            f"{cx:8.3f}{cy + bead_offset:8.3f}{cz:8.3f}"
        )
    lines.append(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}")
    path.write_text("\n".join(lines) + "\n")


def write_top(path: Path, chain_count: int) -> None:
    mol_lines = "\n".join(f"Pctn_{i}    1" for i in range(1, chain_count + 1))
    includes = "\n".join(
        f'#include "toppar_custom/sudowoodo_pectin_{i}.itp"'
        for i in range(1, chain_count + 1)
    )
    text = textwrap.dedent(
        f"""\
        ;;;;;; Bead-count sweep topology

        #include "toppar_custom/sudowoodo_base.itp"
        {includes}

        [ system ]
        Pectin bead-count sweep system

        [ molecules ]
        {mol_lines}
        """
    )
    path.write_text(text)


def write_mdp_files(case_dir: Path, prod_ns: float, dt_ps: float) -> None:
    nsteps = production_steps(prod_ns, dt_ps)
    em = textwrap.dedent(
        """\
        integrator  = steep
        emtol       = 100.0
        emstep      = 0.01
        nsteps      = 50000
        nstlist         = 1
        cutoff-scheme   = Verlet
        ns_type         = grid
        coulombtype     = reaction-field
        vdw-type        = cutoff
        vdw-modifier    = Potential-shift-verlet
        rcoulomb        = 2.5
        rvdw            = 2.5
        pbc             = xyz
        verlet-buffer-tolerance  = 0.005
        """
    )
    eq = textwrap.dedent(
        """\
        integrator               = sd
        dt                       = 0.05
        nsteps                   = 400000
        nstcomm                  = 100
        comm-grps                = system
        nstxout                  = 0
        nstvout                  = 0
        nstfout                  = 0
        nstlog                   = 50000
        nstenergy                = 50000
        nstxout-compressed       = 5000
        compressed-x-precision   = 1000
        compressed-x-grps        = system
        energygrps               = system
        cutoff-scheme            = Verlet
        nstlist                  = 20
        ns_type                  = grid
        pbc                      = xyz
        verlet-buffer-tolerance  = 0.005
        coulombtype              = reaction-field
        rcoulomb                 = 6.0
        epsilon_r                = 15
        epsilon_rf               = 0
        vdw_type                 = cutoff
        vdw-modifier             = Potential-shift-verlet
        rvdw                     = 6.0
        rlist                    = 6.0
        tcoupl                   = v-rescale
        tc-grps                  = system
        tau_t                    = 0.05
        ref_t                    = 300
        Pcoupl                   = no
        gen_vel                  = no
        """
    )
    prod = textwrap.dedent(
        f"""\
        integrator               = sd
        dt                       = {dt_ps}
        nsteps                   = {nsteps}
        nstcomm                  = 100
        nstxout                  = 0
        nstvout                  = 0
        nstfout                  = 0
        nstlog                   = 100000
        nstenergy                = 100000
        nstxout-compressed       = 5000
        compressed-x-precision   = 100
        cutoff-scheme            = Verlet
        nstlist                  = 20
        rlist                    = 20
        rvdw                     = 7
        pbc                      = xyz
        verlet-buffer-tolerance  = 0.005
        coulombtype              = Reaction-Field
        rcoulomb                 = 7
        epsilon_r                = 15
        epsilon_rf               = 0
        vdw_type                 = cutoff
        vdw-modifier             = Potential-shift-verlet
        tcoupl                   = v-rescale
        tc-grps                  = system
        tau_t                    = 0.05
        ref_t                    = 300
        Pcoupl                   = no
        gen_vel                  = no
        gen_temp                 = 300
        """
    )
    (case_dir / "EM.mdp").write_text(em)
    (case_dir / "EQ.mdp").write_text(eq)
    (case_dir / "production.mdp").write_text(prod)


def write_run_sh(case_dir: Path, args: argparse.Namespace) -> None:
    text = textwrap.dedent(
        f"""\
        #!/bin/bash
        set -euo pipefail

        {args.gmx} grompp -f EM.mdp -c afm_system.gro -p afm_system.top -o EM.tpr
        {args.gmx} mdrun -deffnm EM -ddcheck -ntmpi {args.ntmpi} -ntomp {args.ntomp} -dlb no

        {args.gmx} grompp -f EQ.mdp -c EM.gro -p afm_system.top -o EQ.tpr -maxwarn 2
        {args.gmx} mdrun -deffnm EQ -ddcheck -ntmpi {args.ntmpi} -ntomp {args.ntomp} -dlb no -v

        {args.gmx} grompp -f production.mdp -c EQ.gro -p afm_system.top -o production.tpr
        {args.gmx} mdrun -deffnm production -ddcheck -ntmpi {args.ntmpi} -ntomp {args.ntomp} -dlb no -v
        """
    )
    path = case_dir / "run.sh"
    path.write_text(text)
    os.chmod(path, 0o755)


def write_log(path: Path, pr_count: int, pc_count: int, pn_count: int, args: argparse.Namespace) -> None:
    text = textwrap.dedent(
        f"""\
        Bead-count sweep case
        =====================
        PR beads per fiber : {pr_count}
        PC beads per fiber : {pc_count}
        PN beads per fiber : {pn_count}
        Total beads/fiber  : {BEADS_PER_FIBER}
        Chains             : {args.chains}
        PR epsilon (kJ/mol): {args.pr_epsilon}
        PN epsilon (kJ/mol): {args.pn_epsilon}
        PC epsilon (kJ/mol): {args.pc_epsilon}
        Box (nm)           : {args.box[0]} x {args.box[1]} x {args.box[2]}
        Production (ns)    : {args.prod_ns}
        Timestep (ps)      : {args.dt_ps}
        Production steps   : {production_steps(args.prod_ns, args.dt_ps)}
        """
    )
    path.write_text(text)


# ---------------------------------------------------------------------------
# Case builder
# ---------------------------------------------------------------------------

def build_case(
    case_dir: Path,
    pr_count: int,
    pc_count: int,
    epsilon_by_type: dict,
    args: argparse.Namespace,
) -> None:
    if case_dir.exists():
        shutil.rmtree(case_dir)
    toppar_dir = case_dir / "toppar_custom"
    ensure_dir(toppar_dir)

    pn_count = BEADS_PER_FIBER - pr_count - pc_count
    rng = random.Random(args.seed)
    assignments = assign_all_chain_bead_epsilons(
        args.chains,
        rng=rng,
        pc_per_fiber=pc_count,
        pr_per_fiber=pr_count,
        epsilon_by_type=epsilon_by_type,
    )

    chain_positions = generate_chain_positions(args.chains, tuple(args.box), args.spacing)
    write_gro(case_dir / "afm_system.gro", assignments, chain_positions, tuple(args.box))
    write_top(case_dir / "afm_system.top", args.chains)
    write_per_bead_base_itp(toppar_dir / "sudowoodo_base.itp", assignments, epsilon_by_type=epsilon_by_type)
    write_per_fiber_pectin_itps(toppar_dir, assignments)
    write_assignment_report(toppar_dir / "pectin_assignment_report.txt", assignments)
    write_mdp_files(case_dir, args.prod_ns, args.dt_ps)
    write_run_sh(case_dir, args)
    write_log(case_dir / "afm_build.log", pr_count, pc_count, pn_count, args)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    args.box = tuple(args.box)

    epsilon_by_type = pectin_epsilon_by_type(args.pr_epsilon, args.pn_epsilon, args.pc_epsilon)
    ensure_dir(args.out)

    for pr_count, pc_count in iter_bead_counts(args.pr_max, args.pc_max, args.bead_step):
        pn_count = BEADS_PER_FIBER - pr_count - pc_count
        case_name = f"pr{pr_count:02d}_pc{pc_count:02d}_pn{pn_count:02d}"
        case_dir = args.out / case_name
        build_case(case_dir, pr_count, pc_count, epsilon_by_type, args)
        print(f"[ok] prepared {case_dir}  (PR={pr_count}, PC={pc_count}, PN={pn_count})")


if __name__ == "__main__":
    main()
