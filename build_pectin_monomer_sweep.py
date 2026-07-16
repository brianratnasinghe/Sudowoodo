#!/usr/bin/env python3
import argparse
import math
import os
import shutil
import textwrap
from decimal import Decimal
from pathlib import Path

DEFAULT_COUNT = 100
DEFAULT_BOX = (40.0, 40.0, 40.0)
DEFAULT_SPACING = 3.5
DEFAULT_EPSILON_START = Decimal("0.1")
DEFAULT_EPSILON_STOP = Decimal("5.0")
DEFAULT_EPSILON_STEP = Decimal("0.1")
DEFAULT_PRODUCTION_NS = Decimal("100")
DEFAULT_DT_PS = Decimal("0.1")
PECTIN_SIGMA = 1.0
PECTIN_MASS = 26.6


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a monomeric pectin epsilon sweep with 100 unbonded beads."
    )
    parser.add_argument("--out", type=Path, required=True, help="Output folder")
    parser.add_argument("--count", type=int, default=DEFAULT_COUNT, help="Number of pectin beads")
    parser.add_argument("--box", type=float, nargs=3, default=DEFAULT_BOX, metavar=("X", "Y", "Z"))
    parser.add_argument("--spacing", type=float, default=DEFAULT_SPACING, help="Grid spacing between beads")
    parser.add_argument("--epsilon-start", type=Decimal, default=DEFAULT_EPSILON_START)
    parser.add_argument("--epsilon-stop", type=Decimal, default=DEFAULT_EPSILON_STOP)
    parser.add_argument("--epsilon-step", type=Decimal, default=DEFAULT_EPSILON_STEP)
    parser.add_argument("--prod-ns", type=Decimal, default=DEFAULT_PRODUCTION_NS, help="Production length in ns")
    parser.add_argument("--dt-ps", type=Decimal, default=DEFAULT_DT_PS, help="Production timestep in ps")
    parser.add_argument("--gmx", type=str, default="gmx")
    parser.add_argument("--ntomp", type=int, default=24)
    parser.add_argument("--ntmpi", type=int, default=1)
    return parser.parse_args()


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def format_decimal(value):
    return format(value.normalize(), "f")


def iter_epsilon_values(start, stop, step):
    if step <= 0:
        raise ValueError("epsilon-step must be positive")
    if stop < start:
        raise ValueError("epsilon-stop must be greater than or equal to epsilon-start")

    epsilon = start
    values = []
    while epsilon <= stop:
        values.append(epsilon)
        epsilon += step
    return values


def production_steps(prod_ns, dt_ps):
    if prod_ns <= 0:
        raise ValueError("prod-ns must be positive")
    if dt_ps <= 0:
        raise ValueError("dt-ps must be positive")
    total_ps = prod_ns * Decimal("1000")
    return int(total_ps / dt_ps)


def generate_positions(count, box, spacing):
    if count <= 0:
        raise ValueError("count must be positive")
    if spacing <= 0:
        raise ValueError("spacing must be positive")

    nx = math.ceil(math.sqrt(count))
    ny = math.ceil(count / nx)
    x_extent = (nx - 1) * spacing
    y_extent = (ny - 1) * spacing

    if x_extent > box[0] or y_extent > box[1]:
        raise ValueError(
            f"box {box} is too small for {count} beads at spacing {spacing}"
        )

    x0 = (box[0] - x_extent) / 2.0
    y0 = (box[1] - y_extent) / 2.0
    z = box[2] / 2.0

    positions = []
    for bead_idx in range(count):
        ix = bead_idx % nx
        iy = bead_idx // nx
        positions.append((x0 + ix * spacing, y0 + iy * spacing, z))
    return positions


def write_gro(path, count, positions, box):
    lines = ["Monomeric pectin bead system", f"{count:5d}"]
    for atom_id, (x, y, z) in enumerate(positions, start=1):
        lines.append(f"{atom_id:5d}{'Pctn':<5}{'P1':>5}{atom_id:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}")
    Path(path).write_text("\n".join(lines) + "\n")


def write_top(path, count):
    text = textwrap.dedent(
        f"""\
        ;;;;;; Monomeric pectin sweep topology

        #include "toppar_custom/sudowoodo_base.itp"
        #include "toppar_custom/sudowoodo_pectin.itp"

        [ system ]
        Monomeric pectin bead system

        [ molecules ]
        Pctn {count}
        """
    )
    Path(path).write_text(text)


def write_base_itp(path, epsilon):
    eps_text = format_decimal(epsilon)
    text = textwrap.dedent(
        f"""\
        ;;;;; Monomeric pectin force field
        ;;;;; All pectin beads are identical and unbonded in this system.

        [ defaults ]
        1 2 ; sigma-epsilon format of LJ parameters

        [ atomtypes ]
        ; type at.num   mass charge ptype      sigma      epsilon
          P        1   {PECTIN_MASS:4.1f}  0.000     A   {PECTIN_SIGMA:8.6f}   {float(epsilon):10.6f}

        [ nonbond_params ]
        ;   i     j  func  sigma(nm)   epsilon(kJ/mol)
            P     P   1    {PECTIN_SIGMA:8.6f}   {float(epsilon):10.6f}
        ; sweep epsilon_PP = {eps_text} kJ/mol
        """
    )
    Path(path).write_text(text)


def write_pectin_itp(path):
    text = textwrap.dedent(
        """\
        ;;;;; Single-bead pectin monomer

        [moleculetype]
        ; molname nrexcl
          Pctn    1

        [atoms]
        ; id type resnr residu atom cgnr charge
          1  P    1     Pctn   P1   1    0
        """
    )
    Path(path).write_text(text)


def write_mdp_files(case_dir, prod_ns, dt_ps):
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
        dt                       = {format_decimal(dt_ps)}
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


def write_run_sh(case_dir, args):
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


def write_log(path, epsilon, args):
    text = textwrap.dedent(
        f"""\
        Monomeric pectin sweep case
        ===========================
        Pectin bead count: {args.count}
        Pectin-Pectin epsilon (kJ/mol): {format_decimal(epsilon)}
        Box (nm): {args.box[0]} x {args.box[1]} x {args.box[2]}
        Production length (ns): {format_decimal(args.prod_ns)}
        Production timestep (ps): {format_decimal(args.dt_ps)}
        Production steps: {production_steps(args.prod_ns, args.dt_ps)}
        """
    )
    Path(path).write_text(text)


def build_case(case_dir, epsilon, args):
    if case_dir.exists():
        shutil.rmtree(case_dir)
    ensure_dir(case_dir / "toppar_custom")

    positions = generate_positions(args.count, args.box, args.spacing)
    write_gro(case_dir / "afm_system.gro", args.count, positions, args.box)
    write_top(case_dir / "afm_system.top", args.count)
    write_base_itp(case_dir / "toppar_custom" / "sudowoodo_base.itp", epsilon)
    write_pectin_itp(case_dir / "toppar_custom" / "sudowoodo_pectin.itp")
    write_mdp_files(case_dir, args.prod_ns, args.dt_ps)
    write_run_sh(case_dir, args)
    write_log(case_dir / "afm_build.log", epsilon, args)


def main():
    args = parse_args()
    args.box = tuple(args.box)
    ensure_dir(args.out)

    for epsilon in iter_epsilon_values(args.epsilon_start, args.epsilon_stop, args.epsilon_step):
        case_dir = args.out / f"pp_eps_{format_decimal(epsilon)}"
        build_case(case_dir, epsilon, args)
        print(f"[ok] prepared {case_dir}")


if __name__ == "__main__":
    main()
