#!/usr/bin/env python3
from __future__ import annotations

import argparse
import random
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Union

PECTIN_ITP_TEMPLATE = Path(__file__).resolve().parent / "toppar_custom" / "sudowoodo_pectin.itp"

EPSILON_STEP = 0.1
EPSILON_STEP_COUNT = 50  # total steps from 0.1 to 5.0
PECTIN_SIGMA = 1.0

# Fixed composition per pectin fiber
BEADS_PER_FIBER = 30
PC_PER_FIBER = 2
PR_PER_FIBER = 2
PN_PER_FIBER = BEADS_PER_FIBER - PC_PER_FIBER - PR_PER_FIBER  # 26

# Bonded parameters copied into per-fiber ITP files
K_BOND = 317
K_THETA = 8.26
BOND_LENGTH = 3.103

EPSILON_RANGE_BY_TYPE: Dict[str, Tuple[float, float]] = {
    "PR": (0.1, 2.0),
    "PN": (2.1, 4.0),
    "PC": (4.0, 5.0),
}

CORE_ATOMTYPES = (
    ("C",  100.0, 0.0, 0.0),
    ("X",   50.0, 0.0, 0.0),
    ("P",   26.6, 0.0, 0.0),
    ("PN",  26.6, 0.0, 0.0),
    ("PR",  26.6, 0.0, 0.0),
    ("PC",  26.6, 0.0, 0.0),
)

CORE_NONBOND_PARAMS = (
    ("C",  "C",  1, 2.673000, 1.000000),
    ("C",  "X",  1, 2.087000, 1.000000),
    ("C",  "P",  1, 1.837000, 1.000000),
    ("X",  "X",  1, 1.500000, 1.000000),
    ("X",  "P",  1, 1.250000, 1.000000),
    ("P",  "P",  1, 1.000000, 1.000000),
)


Assignment = Dict[str, Union[float, int, str]]
AssignmentMap = Dict[int, Dict[int, Assignment]]


def _epsilon_steps_for_type(bead_type: str) -> List[float]:
    """Return the list of valid epsilon values (step 0.1) for a given bead type."""
    lo, hi = EPSILON_RANGE_BY_TYPE[bead_type]
    n_steps = round((hi - lo) / EPSILON_STEP) + 1
    return [round(lo + i * EPSILON_STEP, 1) for i in range(n_steps)]


def epsilon_step_values() -> List[float]:
    """All valid pectin epsilon values across all types (0.1 to 5.0, step 0.1)."""
    return [round(step * EPSILON_STEP, 1) for step in range(1, EPSILON_STEP_COUNT + 1)]


def epsilon_to_step_index(epsilon: float) -> int:
    lo = EPSILON_RANGE_BY_TYPE["PR"][0]
    hi = EPSILON_RANGE_BY_TYPE["PC"][1]
    if epsilon < lo or epsilon > hi:
        raise ValueError(
            f"Epsilon must be between {lo:.1f} and {hi:.1f}, got {epsilon}"
        )
    scaled = epsilon / EPSILON_STEP
    if abs(scaled - round(scaled)) > 1e-9:
        raise ValueError(f"Pectin epsilon must be a multiple of {EPSILON_STEP:.1f}, got {epsilon}")
    return int(round(scaled))


def step_index_to_epsilon(step_index: int) -> float:
    return round(step_index * EPSILON_STEP, 1)


def classify_pectin_epsilon(epsilon: float) -> str:
    """Classify epsilon into PR/PN/PC.  PC is checked before PN so that 4.0
    (the shared boundary) is always classified as PC."""
    if EPSILON_RANGE_BY_TYPE["PR"][0] <= epsilon <= EPSILON_RANGE_BY_TYPE["PR"][1]:
        return "PR"
    if EPSILON_RANGE_BY_TYPE["PC"][0] <= epsilon <= EPSILON_RANGE_BY_TYPE["PC"][1]:
        return "PC"
    if EPSILON_RANGE_BY_TYPE["PN"][0] <= epsilon <= EPSILON_RANGE_BY_TYPE["PN"][1]:
        return "PN"
    raise ValueError(f"Unsupported pectin epsilon: {epsilon}")


def _epsilon_step_code(epsilon: float) -> str:
    return f"{epsilon_to_step_index(epsilon):02d}"


def _pectin_atomtype_name(chain_index: int, bead_index: int, bead_type: str, epsilon: float) -> str:
    """Build a unique GROMACS atomtype name for a single pectin bead.

    The name encodes the bead class (PR/PN/PC), its epsilon step index,
    the chain it belongs to, and its position within that chain, e.g.
    ``PNe25c3b12``.
    """
    return f"{bead_type}e{_epsilon_step_code(epsilon)}c{chain_index}b{bead_index}"


def assign_all_chain_bead_epsilons(
    chain_count: int,
    rng: random.Random | None = None,
) -> AssignmentMap:
    """Assign per-bead epsilons for *chain_count* pectin chains.

    Each chain receives a fixed composition of
    ``PC_PER_FIBER`` PC beads, ``PR_PER_FIBER`` PR beads, and
    ``PN_PER_FIBER`` PN beads, placed in a randomised order.
    The epsilon for each bead is drawn uniformly from the 0.1-step
    grid within its type's valid range.
    """
    rng = rng or random.Random()
    assignments: AssignmentMap = {}
    for chain_index in range(1, chain_count + 1):
        bead_types: List[str] = (
            ["PC"] * PC_PER_FIBER + ["PR"] * PR_PER_FIBER + ["PN"] * PN_PER_FIBER
        )
        rng.shuffle(bead_types)
        chain_assignments: Dict[int, Assignment] = {}
        for bead_index, bead_type in enumerate(bead_types, 1):
            epsilon = rng.choice(_epsilon_steps_for_type(bead_type))
            chain_assignments[bead_index] = {
                "chain_index": chain_index,
                "bead_index": bead_index,
                "epsilon": epsilon,
                "bead_type": bead_type,
                "atomtype": _pectin_atomtype_name(chain_index, bead_index, bead_type, epsilon),
            }
        assignments[chain_index] = chain_assignments
    return assignments


def iter_assignments(assignments: AssignmentMap) -> Iterable[Assignment]:
    for chain_index in sorted(assignments):
        for bead_index in sorted(assignments[chain_index]):
            yield assignments[chain_index][bead_index]


def sorted_assignments_by_epsilon(assignments: AssignmentMap) -> List[Assignment]:
    return sorted(
        iter_assignments(assignments),
        key=lambda item: (item["epsilon"], item["chain_index"], item["bead_index"]),
    )


def _base_header() -> str:
    return """;;;;; Sudowoodo Force Field: Particle definition and LJ interactions
;;;;;
;;;;; File generated by build_sweep.py
;;;;;
;;;;; Version: 0.0.1
;;;;;

[ defaults ]
1 2 ; sigma-epsilon format of LJ parameters

"""


def write_per_bead_base_itp(output_path: Path, assignments: AssignmentMap) -> None:
    """Write a complete sudowoodo_base.itp with core types AND per-bead pectin atomtypes."""
    lines = [_base_header(), "[ atomtypes ]", "; name  at.num  mass  charge  ptype  sigma  epsilon"]
    for name, mass, sigma, epsilon in CORE_ATOMTYPES:
        lines.append(f"{name:<12} 1 {mass:>7.1f} 0.000 A {sigma:>9.6f} {epsilon:>11.6f}")
    # One atomtype entry per bead across all chains
    for assignment in iter_assignments(assignments):
        atype = str(assignment["atomtype"])
        eps = float(assignment["epsilon"])
        lines.append(f"{atype:<24} 1   26.6 0.000 A  {PECTIN_SIGMA:.6f} {eps:11.6f}")

    lines.extend([
        "",
        "[ nonbond_params ]",
        ";   i     j  func  sigma(nm)   epsilon(kJ/mol)",
    ])
    for left, right, func, sigma, epsilon in CORE_NONBOND_PARAMS:
        lines.append(f"{left:>12} {right:>12} {func:>3} {sigma:>12.6f} {epsilon:>12.6f}")

    output_path.write_text("\n".join(lines) + "\n")


def append_per_bead_atomtypes(base_itp_path: Path, assignments: AssignmentMap) -> None:
    """Append per-bead pectin atomtypes to an *existing* sudowoodo_base.itp.

    This is used when ``afm_build_sweep.py`` has already written a base ITP
    with custom C/X/P epsilon values; we inject the per-bead lines at the end
    of the ``[ atomtypes ]`` block without disturbing the rest of the file.
    """
    new_entries = []
    for assignment in iter_assignments(assignments):
        atype = str(assignment["atomtype"])
        eps = float(assignment["epsilon"])
        new_entries.append(f"{atype:<24} 1   26.6 0.000 A  {PECTIN_SIGMA:.6f} {eps:11.6f}")

    lines = base_itp_path.read_text().splitlines()
    # Insert the new atomtype lines just before the [ nonbond_params ] section
    insert_at = None
    for i, line in enumerate(lines):
        if re.match(r"^\s*\[\s*nonbond_params", line):
            insert_at = i
            break

    if insert_at is not None:
        lines = lines[:insert_at] + new_entries + [""] + lines[insert_at:]
    else:
        lines = lines + new_entries

    base_itp_path.write_text("\n".join(lines) + "\n")


def write_per_fiber_pectin_itp(
    output_path: Path,
    chain_index: int,
    chain_assignments: Dict[int, Assignment],
) -> None:
    """Generate a per-fiber sudowoodo_pectin_N.itp.

    The molecule is named ``Pctn_<chain_index>``.  Each of the 30 atoms gets
    a unique atomtype name that encodes its class, epsilon, chain and bead
    position.  Bond and angle topology mirrors the template ITP.
    """
    mol_name = f"Pctn_{chain_index}"
    n = len(chain_assignments)
    lines = [
        f"#define k_bond {K_BOND}",
        f"#define k_theta {K_THETA}",
        "",
        "[ moleculetype ]",
        "; molname      nrexcl",
        f"  {mol_name:<14} 1",
        "",
        "[ atoms ]",
        "; id    type                      resnr  residu  atom   cgnr  charge",
    ]
    for bead_index in range(1, n + 1):
        atype = str(chain_assignments[bead_index]["atomtype"])
        lines.append(
            f"  {bead_index:<5} {atype:<24} 1      Pctn    P{bead_index:<5} {bead_index:<5} 0"
        )
    lines += [
        "",
        "[ bonds ]",
        ";  i    j      funct   length  force.c.",
    ]
    for i in range(1, n):
        lines.append(f"  {i}    {i + 1}    1       {BOND_LENGTH}    k_bond")
    lines += [
        "",
        "[ angles ]",
        ";  i    j    k      funct   angle   force.c.",
    ]
    for i in range(1, n - 1):
        lines.append(f"  {i}    {i + 1}    {i + 2}    2       180.0   k_theta")
    output_path.write_text("\n".join(lines) + "\n")


def write_per_fiber_pectin_itps(toppar_dir: Path, assignments: AssignmentMap) -> None:
    """Write one sudowoodo_pectin_N.itp per chain into *toppar_dir*."""
    for chain_index, chain_assignments in assignments.items():
        write_per_fiber_pectin_itp(
            toppar_dir / f"sudowoodo_pectin_{chain_index}.itp",
            chain_index,
            chain_assignments,
        )


def write_assignment_report(output_path: Path, assignments: AssignmentMap) -> None:
    lines = []
    for assignment in sorted_assignments_by_epsilon(assignments):
        lines.append(
            f"{assignment['atomtype']:<24} 1 26.6 0.000 A {PECTIN_SIGMA:.6f} {assignment['epsilon']:.6f}"
        )
    output_path.write_text("\n".join(lines) + "\n")


def build_variant(
    output_dir: Path,
    chain_count: int = 1,
    rng: random.Random | None = None,
) -> AssignmentMap:
    """Build a complete pectin variant in *output_dir*.

    Writes:

    * ``sudowoodo_base.itp`` – core atomtypes + per-bead atomtypes + nonbond_params
    * ``sudowoodo_pectin_N.itp`` for each chain (N = 1 … chain_count)
    * ``pectin_assignment_report.txt`` – sorted assignment listing
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    assignments = assign_all_chain_bead_epsilons(chain_count, rng=rng)
    write_per_bead_base_itp(output_dir / "sudowoodo_base.itp", assignments)
    write_per_fiber_pectin_itps(output_dir, assignments)
    write_assignment_report(output_dir / "pectin_assignment_report.txt", assignments)
    return assignments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate per-fiber pectin ITP files with randomised PR/PN/PC beads.")
    parser.add_argument("--out", type=Path, required=True, help="Output directory")
    parser.add_argument("--chains", type=int, default=1, help="Number of pectin chains (default 1)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed) if args.seed is not None else random.Random()
    build_variant(args.out, chain_count=args.chains, rng=rng)


if __name__ == "__main__":
    main()
