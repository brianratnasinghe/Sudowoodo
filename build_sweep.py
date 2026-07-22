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
# Note: PN and PC both include ε=4.0 in their ranges.  classify_pectin_epsilon()
# resolves the overlap by checking PC before PN, so ε=4.0 is classified as PC.
# In the catalog, "PN2.0" (absolute ε=4.0, PN offset 2.0) and "PC4.0" are
# distinct named types; each can appear independently in a pectin fiber.

CORE_ATOMTYPES = (
    ("C",  100.0, 0.0, 0.0),
    ("X",   50.0, 0.0, 0.0),
    ("P",   26.6, 0.0, 0.0),
)

CORE_NONBOND_PARAMS = (
    ("C",  "C",  1, 2.673000, 1.000000),
    ("C",  "X",  1, 2.087000, 1.000000),
    ("C",  "P",  1, 1.837000, 1.000000),
    ("X",  "X",  1, 1.500000, 1.000000),
    ("X",  "P",  1, 1.250000, 1.000000),
    ("P",  "P",  1, 1.000000, 1.000000),
)

# Sigma values for core beads interacting with catalog pectin beads.
_CORE_CATALOG_SIGMA: Tuple[Tuple[str, float], ...] = (
    ("C", 1.837),
    ("X", 1.250),
    ("P", 1.000),
)
# Epsilon for each core bead type interacting with catalog pectin beads.
# C and X default to 2.5; P remains at 1.0 so PP can be treated as a fallback
# in AFM remapping without overwriting assigned pectin catalog interactions.
_CORE_CATALOG_EPSILON: Dict[str, float] = {"C": 2.5, "X": 2.5, "P": 1.0}


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


# ---------------------------------------------------------------------------
# Catalog type naming
# ---------------------------------------------------------------------------

def _catalog_type_name(bead_type: str, epsilon: float) -> str:
    """Return the shared catalog atomtype name for a bead class + epsilon.

    Each bead type has a descriptive prefix and a zero-padded 1-based step
    index counted from the low end of the type's epsilon range:

    * ``PctRep`` (Pectin Repulsive)  steps 01–20 for ε 0.1–2.0
    * ``PctNeu`` (Pectin Neutral)    steps 01–20 for ε 2.1–4.0
    * ``PctXlk`` (Pectin Crosslink)  steps 01–11 for ε 4.0–5.0
    """
    lo = EPSILON_RANGE_BY_TYPE[bead_type][0]
    step = round((epsilon - lo) / EPSILON_STEP) + 1
    if bead_type == "PR":
        return f"PctRep{step:02d}"
    if bead_type == "PN":
        return f"PctNeu{step:02d}"
    if bead_type == "PC":
        return f"PctXlk{step:02d}"
    raise ValueError(f"Unknown bead type: {bead_type!r}")


def catalog_type_names() -> List[str]:
    """Return all catalog pectin type names in PR → PN → PC order."""
    names: List[str] = []
    for btype in ("PR", "PN", "PC"):
        for eps in _epsilon_steps_for_type(btype):
            names.append(_catalog_type_name(btype, eps))
    return names


def _catalog_items() -> List[Tuple[str, str, float]]:
    """Return ``(bead_type, type_name, abs_epsilon)`` for every catalog type."""
    items: List[Tuple[str, str, float]] = []
    for btype in ("PR", "PN", "PC"):
        for eps in _epsilon_steps_for_type(btype):
            items.append((btype, _catalog_type_name(btype, eps), eps))
    return items


def _catalog_stored_eps(bead_type: str, epsilon: float) -> float:
    """Interaction epsilon stored in ``[nonbond_params]`` for a catalog type.

    * PR: absolute epsilon (same as offset from 0)
    * PN: epsilon − 2.0 (offset from PN base)
    * PC: absolute epsilon
    """
    if bead_type == "PN":
        return round(epsilon - 2.0, 1)
    return round(epsilon, 1)


def _catalog_cross_eps(btype_i: str, eps_i: float, btype_j: str, eps_j: float) -> float:
    """Compute the ``[nonbond_params]`` epsilon for a pair of catalog types.

    Combining rules (matching the problem-statement examples):

    * Same-class **PC** pairs: arithmetic mean of absolute epsilons
      (e.g. PC4.0 × PC4.1 → 4.05).
    * Same-class **PN** pairs: max of stored offsets (epsilon − 2.0)
      (e.g. PN.1 × PN.2 → 0.2 = max(0.1, 0.2)).
    * Same-class **PR** pairs: max of stored absolute epsilons
      (e.g. PR.1 × PR.2 → 0.2 = max(0.1, 0.2)).
    * Cross-class pairs (PR×PN, PR×PC, PN×PC): arithmetic mean of stored
      epsilons.  Note that PR and PN stored values use different scales (PR
      absolute, PN offset from 2.0), so the cross-class result is a
      scale-mixed average; this is acceptable because the problem statement
      does not prescribe cross-class interactions explicitly.
    """
    stored_i = _catalog_stored_eps(btype_i, eps_i)
    stored_j = _catalog_stored_eps(btype_j, eps_j)
    if btype_i == btype_j == "PC":
        return round((stored_i + stored_j) / 2.0, 6)
    if btype_i == btype_j:  # PR×PR or PN×PN
        return round(max(stored_i, stored_j), 6)
    # Cross-class: arithmetic mean of stored epsilons
    return round((stored_i + stored_j) / 2.0, 6)


def _pectin_atomtype_name(chain_index: int, bead_index: int, bead_type: str, epsilon: float) -> str:
    """Return the shared catalog atomtype name for a pectin bead.

    ``chain_index`` and ``bead_index`` are accepted for API compatibility but
    are not encoded in the name; multiple beads may share the same catalog type.
    """
    return _catalog_type_name(bead_type, epsilon)


def assign_all_chain_bead_epsilons(
    chain_count: int,
    rng: random.Random | None = None,
    pc_per_fiber: int = PC_PER_FIBER,
    pr_per_fiber: int = PR_PER_FIBER,
) -> AssignmentMap:
    """Assign per-bead epsilons for *chain_count* pectin chains.

    Each chain receives a composition of *pc_per_fiber* PC beads,
    *pr_per_fiber* PR beads, and the remainder as PN beads (up to
    ``BEADS_PER_FIBER`` total), placed in a randomised order.
    The epsilon for each bead is drawn uniformly from the 0.1-step
    grid within its type's valid range.
    """
    pn = BEADS_PER_FIBER - pc_per_fiber - pr_per_fiber
    if pn < 0:
        raise ValueError(
            f"pc_per_fiber ({pc_per_fiber}) + pr_per_fiber ({pr_per_fiber}) "
            f"exceeds BEADS_PER_FIBER ({BEADS_PER_FIBER})"
        )
    rng = rng or random.Random()
    assignments: AssignmentMap = {}
    for chain_index in range(1, chain_count + 1):
        bead_types: List[str] = (
            ["PC"] * pc_per_fiber + ["PR"] * pr_per_fiber + ["PN"] * pn
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
    """Write a complete ``sudowoodo_base.itp`` with core types, the full catalog
    of pectin atomtypes, and all pairwise ``[nonbond_params]``.

    The *assignments* parameter is accepted for API compatibility; the catalog
    written is always the complete set of PR/PN/PC types regardless of which
    specific types appear in the assignments.
    """
    lines = [_base_header(), "[ atomtypes ]", "; name  at.num  mass  charge  ptype  sigma  epsilon"]
    for name, mass, sigma, epsilon in CORE_ATOMTYPES:
        lines.append(f"{name:<12} 1 {mass:>7.1f} 0.000 A {sigma:>9.6f} {epsilon:>11.6f}")
    # Catalog pectin types — interactions defined entirely via nonbond_params
    for _btype, cname, _eps in _catalog_items():
        lines.append(f"{cname:<16} 1   26.6 0.000 A   0.000000    0.000000")

    lines.extend([
        "",
        "[ nonbond_params ]",
        ";   i     j  func  sigma(nm)   epsilon(kJ/mol)",
    ])
    # Core C/X/P pairs (unchanged)
    for left, right, func, sigma, epsilon in CORE_NONBOND_PARAMS:
        lines.append(f"{left:>12} {right:>12} {func:>3} {sigma:>12.6f} {epsilon:>12.6f}")
    # Core beads × catalog pectin types
    items = _catalog_items()
    for core, core_sigma in _CORE_CATALOG_SIGMA:
        for _btype, cname, _eps in items:
            lines.append(f"{core:>12} {cname:>16} {1:>3} {core_sigma:>12.6f} {_CORE_CATALOG_EPSILON[core]:>12.6f}")
    # All catalog × catalog pairs (upper triangle, i ≤ j)
    n = len(items)
    for i in range(n):
        btype_i, name_i, eps_i = items[i]
        for j in range(i, n):
            btype_j, name_j, eps_j = items[j]
            cross = _catalog_cross_eps(btype_i, eps_i, btype_j, eps_j)
            lines.append(f"{name_i:>16} {name_j:>16} {1:>3} {PECTIN_SIGMA:>12.6f} {cross:>12.6f}")

    output_path.write_text("\n".join(lines) + "\n")


def append_per_bead_atomtypes(base_itp_path: Path, assignments: AssignmentMap) -> None:
    """Inject the full catalog of pectin atomtypes into an *existing* base ITP.

    When ``afm_build_sweep.py`` writes a custom base ITP (with bespoke C/X/P
    parameters), this function inserts the complete PR/PN/PC catalog just before
    the ``[nonbond_params]`` block so that the per-fiber pectin ITPs can
    reference the shared catalog types.

    If all catalog types are already present in the file the function is a
    no-op to avoid duplication.
    """
    existing_text = base_itp_path.read_text()
    # Build list of catalog entries not yet present
    new_entries = []
    for _btype, cname, _eps in _catalog_items():
        if cname not in existing_text:
            new_entries.append(f"{cname:<16} 1   26.6 0.000 A   0.000000    0.000000")

    if not new_entries:
        return  # catalog already fully present

    lines = existing_text.splitlines()
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

    The molecule is named ``Pctn_<chain_index>``.  Each of the 30 atoms
    references a shared catalog atomtype from the base ITP.  Bond and angle
    topology mirrors the template ITP.
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


def write_default_pectin_itp(output_path: Path, seed: int = 0) -> None:
    """Write a static ``sudowoodo_pectin.itp`` template using catalog atomtypes.

    The molecule is named ``Pctn`` and uses the standard composition of
    ``PC_PER_FIBER`` PctXlk + ``PR_PER_FIBER`` PctRep + ``PN_PER_FIBER`` PctNeu
    beads drawn with the given *seed* for reproducibility.  All atomtype names
    reference shared catalog entries defined in ``sudowoodo_base.itp``.
    """
    assignments = assign_all_chain_bead_epsilons(1, rng=random.Random(seed))
    chain_assignments = assignments[1]
    n = len(chain_assignments)
    lines = [
        f"#define k_bond {K_BOND}",
        f"#define k_theta {K_THETA}",
        "",
        "[ moleculetype ]",
        "; molname      nrexcl",
        "  Pctn         1",
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
    pc_per_fiber: int = PC_PER_FIBER,
    pr_per_fiber: int = PR_PER_FIBER,
) -> AssignmentMap:
    """Build a complete pectin variant in *output_dir*.

    Writes:

    * ``sudowoodo_base.itp`` – core C/X/P atomtypes + full PR/PN/PC catalog +
      all pairwise nonbond_params
    * ``sudowoodo_pectin_N.itp`` for each chain (N = 1 … chain_count)
    * ``pectin_assignment_report.txt`` – sorted assignment listing
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    assignments = assign_all_chain_bead_epsilons(chain_count, rng=rng,
                                                 pc_per_fiber=pc_per_fiber,
                                                 pr_per_fiber=pr_per_fiber)
    write_per_bead_base_itp(output_dir / "sudowoodo_base.itp", assignments)
    write_per_fiber_pectin_itps(output_dir, assignments)
    write_assignment_report(output_dir / "pectin_assignment_report.txt", assignments)
    return assignments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate per-fiber pectin ITP files with randomised PR/PN/PC beads.")
    parser.add_argument("--out", type=Path, required=True, help="Output directory")
    parser.add_argument("--chains", type=int, default=1, help="Number of pectin chains (default 1)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--pc", type=int, default=PC_PER_FIBER,
                        help=f"Crosslinking (PC) beads per fiber (default {PC_PER_FIBER})")
    parser.add_argument("--pr", type=int, default=PR_PER_FIBER,
                        help=f"Repulsion (PR) beads per fiber (default {PR_PER_FIBER})")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed) if args.seed is not None else random.Random()
    build_variant(args.out, chain_count=args.chains, rng=rng,
                  pc_per_fiber=args.pc, pr_per_fiber=args.pr)


if __name__ == "__main__":
    main()
