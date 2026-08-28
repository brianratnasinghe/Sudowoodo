#!/usr/bin/env python3
from __future__ import annotations

import argparse
import random
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Union

PECTIN_ITP_TEMPLATE = Path(__file__).resolve().parent / "toppar_custom" / "sudowoodo_pectin.itp"

PECTIN_SIGMA = 1.0

BEADS_PER_FIBER = 30
PC_PER_FIBER = 2
PR_PER_FIBER = 2
PN_PER_FIBER = BEADS_PER_FIBER - PC_PER_FIBER - PR_PER_FIBER

DEFAULT_EPSILON_BY_TYPE: Dict[str, float] = {
    "PR": 0.4,
    "PN": 2.2,
    "PC": 4.8,
}

DEFAULT_PCT_CROSS_EPSILON: float = 2.5

PECTIN_TYPE_NAMES: Dict[str, str] = {
    "PR": "PctRep",
    "PN": "PctNeu",
    "PC": "PctXlk",
}

K_BOND = 317
K_THETA = 8.26
BOND_LENGTH = 3.103

CORE_ATOMTYPES = (
    ("C", 100.0, 0.0, 0.0),
    ("X", 50.0, 0.0, 0.0),
    ("P", 26.6, 0.0, 0.0),
)

CORE_NONBOND_PARAMS = (
    ("C", "C", 1, 2.673000, 2.500000),
    ("C", "X", 1, 2.087000, 2.500000),
    ("C", "P", 1, 1.837000, 2.500000),
    ("X", "X", 1, 1.500000, 2.500000),
    ("X", "P", 1, 1.250000, 2.500000),
    ("P", "P", 1, 1.000000, 2.500000),
)

CORE_CATALOG_SIGMA: Tuple[Tuple[str, float], ...] = (
    ("C", 1.837),
    ("X", 1.250),
    ("P", 1.000),
)

CORE_CATALOG_EPSILON: Dict[str, float] = {"C": 2.5, "X": 2.5, "P": 2.5}

Assignment = Dict[str, Union[float, int, str]]
AssignmentMap = Dict[int, Dict[int, Assignment]]


def pectin_epsilon_by_type(
    pr_epsilon: float = DEFAULT_EPSILON_BY_TYPE["PR"],
    pn_epsilon: float = DEFAULT_EPSILON_BY_TYPE["PN"],
    pc_epsilon: float = DEFAULT_EPSILON_BY_TYPE["PC"],
) -> Dict[str, float]:
    values = {
        "PR": float(pr_epsilon),
        "PN": float(pn_epsilon),
        "PC": float(pc_epsilon),
    }
    for bead_type, epsilon in values.items():
        if epsilon < 0:
            raise ValueError(f"{bead_type} epsilon must be non-negative, got {epsilon}")
    return values


def _resolve_pectin_epsilon_by_type(
    epsilon_by_type: Dict[str, float] | None = None,
    pr_epsilon: float | None = None,
    pn_epsilon: float | None = None,
    pc_epsilon: float | None = None,
) -> Dict[str, float]:
    values = dict(DEFAULT_EPSILON_BY_TYPE)
    if epsilon_by_type:
        values.update({key: float(value) for key, value in epsilon_by_type.items() if key in values})
    if pr_epsilon is not None:
        values["PR"] = float(pr_epsilon)
    if pn_epsilon is not None:
        values["PN"] = float(pn_epsilon)
    if pc_epsilon is not None:
        values["PC"] = float(pc_epsilon)
    return pectin_epsilon_by_type(values["PR"], values["PN"], values["PC"])


def catalog_type_names() -> List[str]:
    return [PECTIN_TYPE_NAMES[bead_type] for bead_type in ("PR", "PN", "PC")]


def is_pectin_atomtype(name: str) -> bool:
    return name.startswith("Pct")


def _catalog_items(epsilon_by_type: Dict[str, float]) -> List[Tuple[str, str, float]]:
    return [
        (bead_type, PECTIN_TYPE_NAMES[bead_type], float(epsilon_by_type[bead_type]))
        for bead_type in ("PR", "PN", "PC")
    ]


def _catalog_cross_eps(btype_i: str, eps_i: float, btype_j: str, eps_j: float, cross_epsilon: float | None = None) -> float:
    if btype_i == btype_j:
        return round(float(eps_i), 6)
    if cross_epsilon is not None:
        return round(float(cross_epsilon), 6)
    return round((float(eps_i) + float(eps_j)) / 2.0, 6)


def _pectin_atomtype_name(chain_index: int, bead_index: int, bead_type: str, epsilon: float) -> str:
    del chain_index, bead_index, epsilon
    return PECTIN_TYPE_NAMES[bead_type]


def assign_all_chain_bead_epsilons(
    chain_count: int,
    rng: random.Random | None = None,
    pc_per_fiber: int = PC_PER_FIBER,
    pr_per_fiber: int = PR_PER_FIBER,
    epsilon_by_type: Dict[str, float] | None = None,
) -> AssignmentMap:
    pn_per_fiber = BEADS_PER_FIBER - pc_per_fiber - pr_per_fiber
    if pn_per_fiber < 0:
        raise ValueError(
            f"pc_per_fiber ({pc_per_fiber}) + pr_per_fiber ({pr_per_fiber}) exceeds "
            f"BEADS_PER_FIBER ({BEADS_PER_FIBER})"
        )

    epsilon_by_type = _resolve_pectin_epsilon_by_type(epsilon_by_type)
    rng = rng or random.Random()
    assignments: AssignmentMap = {}

    for chain_index in range(1, chain_count + 1):
        bead_types: List[str] = ["PC"] * pc_per_fiber + ["PR"] * pr_per_fiber + ["PN"] * pn_per_fiber
        rng.shuffle(bead_types)
        chain_assignments: Dict[int, Assignment] = {}
        for bead_index, bead_type in enumerate(bead_types, 1):
            epsilon = epsilon_by_type[bead_type]
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


def write_per_bead_base_itp(
    output_path: Path,
    assignments: AssignmentMap,
    epsilon_by_type: Dict[str, float] | None = None,
    core_nonbond_params: Tuple[Tuple[str, str, int, float, float], ...] = CORE_NONBOND_PARAMS,
    core_catalog_epsilon: Dict[str, float] | None = None,
    pct_cross_epsilon: float | None = DEFAULT_PCT_CROSS_EPSILON,
) -> None:
    del assignments
    epsilon_by_type = _resolve_pectin_epsilon_by_type(epsilon_by_type)
    core_catalog_epsilon = dict(CORE_CATALOG_EPSILON if core_catalog_epsilon is None else core_catalog_epsilon)
    items = _catalog_items(epsilon_by_type)

    lines = [_base_header(), "[ atomtypes ]", "; name  at.num  mass  charge  ptype  sigma  epsilon"]
    for name, mass, sigma, epsilon in CORE_ATOMTYPES:
        lines.append(f"{name:<12} 1 {mass:>7.1f} 0.000 A {sigma:>9.6f} {epsilon:>11.6f}")
    for _bead_type, atomtype, _epsilon in items:
        lines.append(f"{atomtype:<16} 1   26.6 0.000 A   0.000000    0.000000")

    lines.extend([
        "",
        "[ nonbond_params ]",
        ";   i     j  func  sigma(nm)   epsilon(kJ/mol)",
    ])
    for left, right, func, sigma, epsilon in core_nonbond_params:
        lines.append(f"{left:>12} {right:>12} {func:>3} {sigma:>12.6f} {epsilon:>12.6f}")
    for core, core_sigma in CORE_CATALOG_SIGMA:
        for _bead_type, atomtype, _epsilon in items:
            lines.append(
                f"{core:>12} {atomtype:>16} {1:>3} {core_sigma:>12.6f} {float(core_catalog_epsilon[core]):>12.6f}"
            )
    for i, (btype_i, name_i, eps_i) in enumerate(items):
        for btype_j, name_j, eps_j in items[i:]:
            cross = _catalog_cross_eps(btype_i, eps_i, btype_j, eps_j, pct_cross_epsilon)
            lines.append(f"{name_i:>16} {name_j:>16} {1:>3} {PECTIN_SIGMA:>12.6f} {cross:>12.6f}")

    output_path.write_text("\n".join(lines) + "\n")


def append_per_bead_atomtypes(
    base_itp_path: Path,
    assignments: AssignmentMap,
    epsilon_by_type: Dict[str, float] | None = None,
    pct_cross_epsilon: float | None = DEFAULT_PCT_CROSS_EPSILON,
) -> None:
    text = base_itp_path.read_text() if base_itp_path.exists() else ""
    core_nonbond_params = CORE_NONBOND_PARAMS
    core_catalog_epsilon = CORE_CATALOG_EPSILON
    if text:
        core_nonbond_params, core_catalog_epsilon = parse_existing_base_itp_settings(base_itp_path)
    write_per_bead_base_itp(
        base_itp_path,
        assignments,
        epsilon_by_type=epsilon_by_type,
        core_nonbond_params=core_nonbond_params,
        core_catalog_epsilon=core_catalog_epsilon,
        pct_cross_epsilon=pct_cross_epsilon,
    )


def parse_existing_base_itp_settings(
    base_itp_path: Path,
) -> Tuple[Tuple[Tuple[str, str, int, float, float], ...], Dict[str, float]]:
    text = base_itp_path.read_text()
    core_pairs: Dict[Tuple[str, str], Tuple[int, float, float]] = {}
    core_catalog_epsilon = dict(CORE_CATALOG_EPSILON)
    in_nonbond = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("["):
            in_nonbond = line == "[ nonbond_params ]"
            continue
        if not in_nonbond:
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        left, right = parts[0], parts[1]
        func = int(parts[2])
        sigma = float(parts[3])
        epsilon = float(parts[4])
        if left in {"C", "X", "P"} and right in {"C", "X", "P"}:
            core_pairs[(left, right)] = (func, sigma, epsilon)
        elif left in {"C", "X", "P"} and is_pectin_atomtype(right):
            core_catalog_epsilon[left] = epsilon
        elif right in {"C", "X", "P"} and is_pectin_atomtype(left):
            core_catalog_epsilon[right] = epsilon

    parsed_core_nonbond_params = []
    for left, right, default_func, default_sigma, default_epsilon in CORE_NONBOND_PARAMS:
        func, sigma, epsilon = core_pairs.get((left, right), (default_func, default_sigma, default_epsilon))
        parsed_core_nonbond_params.append((left, right, func, sigma, epsilon))
    return tuple(parsed_core_nonbond_params), core_catalog_epsilon


def _write_pectin_itp_lines(mol_name: str, chain_assignments: Dict[int, Assignment]) -> List[str]:
    n_beads = len(chain_assignments)
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
    for bead_index in range(1, n_beads + 1):
        atomtype = str(chain_assignments[bead_index]["atomtype"])
        bead_type = str(chain_assignments[bead_index]["bead_type"])
        lines.append(f"  {bead_index:<5} {atomtype:<24} 1      Pctn    {bead_type:<6} {bead_index:<5} 0")
    lines.extend([
        "",
        "[ bonds ]",
        ";  i    j      funct   length  force.c.",
    ])
    for bead_index in range(1, n_beads):
        lines.append(f"  {bead_index}    {bead_index + 1}    1       {BOND_LENGTH}    k_bond")
    lines.extend([
        "",
        "[ angles ]",
        ";  i    j    k      funct   angle   force.c.",
    ])
    for bead_index in range(1, n_beads - 1):
        lines.append(f"  {bead_index}    {bead_index + 1}    {bead_index + 2}    2       180.0   k_theta")
    return lines


def write_per_fiber_pectin_itp(output_path: Path, chain_index: int, chain_assignments: Dict[int, Assignment]) -> None:
    output_path.write_text("\n".join(_write_pectin_itp_lines(f"Pctn_{chain_index}", chain_assignments)) + "\n")


def write_per_fiber_pectin_itps(toppar_dir: Path, assignments: AssignmentMap) -> None:
    for chain_index, chain_assignments in assignments.items():
        write_per_fiber_pectin_itp(toppar_dir / f"sudowoodo_pectin_{chain_index}.itp", chain_index, chain_assignments)


def write_default_pectin_itp(
    output_path: Path,
    seed: int = 0,
    pc_per_fiber: int = PC_PER_FIBER,
    pr_per_fiber: int = PR_PER_FIBER,
    epsilon_by_type: Dict[str, float] | None = None,
) -> None:
    assignments = assign_all_chain_bead_epsilons(
        1,
        rng=random.Random(seed),
        pc_per_fiber=pc_per_fiber,
        pr_per_fiber=pr_per_fiber,
        epsilon_by_type=epsilon_by_type,
    )
    output_path.write_text("\n".join(_write_pectin_itp_lines("Pctn", assignments[1])) + "\n")


def write_assignment_report(output_path: Path, assignments: AssignmentMap) -> None:
    lines = []
    for assignment in iter_assignments(assignments):
        lines.append(
            f"chain={assignment['chain_index']} bead={assignment['bead_index']} "
            f"type={assignment['bead_type']} atomtype={assignment['atomtype']} epsilon={float(assignment['epsilon']):.6f}"
        )
    output_path.write_text("\n".join(lines) + "\n")


def build_variant(
    output_dir: Path,
    chain_count: int = 1,
    rng: random.Random | None = None,
    pc_per_fiber: int = PC_PER_FIBER,
    pr_per_fiber: int = PR_PER_FIBER,
    epsilon_by_type: Dict[str, float] | None = None,
) -> AssignmentMap:
    output_dir.mkdir(parents=True, exist_ok=True)
    assignments = assign_all_chain_bead_epsilons(
        chain_count,
        rng=rng,
        pc_per_fiber=pc_per_fiber,
        pr_per_fiber=pr_per_fiber,
        epsilon_by_type=epsilon_by_type,
    )
    write_per_bead_base_itp(output_dir / "sudowoodo_base.itp", assignments, epsilon_by_type=epsilon_by_type)
    write_per_fiber_pectin_itps(output_dir, assignments)
    write_assignment_report(output_dir / "pectin_assignment_report.txt", assignments)
    return assignments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate per-fiber pectin ITP files with randomized PR/PN/PC beads.")
    parser.add_argument("--out", type=Path, required=True, help="Output directory")
    parser.add_argument("--chains", type=int, default=1, help="Number of pectin chains (default 1)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--pc", type=int, default=PC_PER_FIBER, help=f"Crosslinking (PC) beads per fiber (default {PC_PER_FIBER})")
    parser.add_argument("--pr", type=int, default=PR_PER_FIBER, help=f"Repulsion (PR) beads per fiber (default {PR_PER_FIBER})")
    parser.add_argument("--pr-epsilon", type=float, default=DEFAULT_EPSILON_BY_TYPE["PR"], help=f"Repulsive bead epsilon (default {DEFAULT_EPSILON_BY_TYPE['PR']})")
    parser.add_argument("--pn-epsilon", type=float, default=DEFAULT_EPSILON_BY_TYPE["PN"], help=f"Neutral bead epsilon (default {DEFAULT_EPSILON_BY_TYPE['PN']})")
    parser.add_argument("--pc-epsilon", type=float, default=DEFAULT_EPSILON_BY_TYPE["PC"], help=f"Crosslink bead epsilon (default {DEFAULT_EPSILON_BY_TYPE['PC']})")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed) if args.seed is not None else random.Random()
    build_variant(
        args.out,
        chain_count=args.chains,
        rng=rng,
        pc_per_fiber=args.pc,
        pr_per_fiber=args.pr,
        epsilon_by_type=pectin_epsilon_by_type(args.pr_epsilon, args.pn_epsilon, args.pc_epsilon),
    )


if __name__ == "__main__":
    main()
