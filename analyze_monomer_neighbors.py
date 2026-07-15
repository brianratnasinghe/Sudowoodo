#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT_DIR = Path(".")
CASE_GLOB = "pp_eps_*"
TOPOLOGY_FILENAME = "afm_system.gro"
TRAJECTORY_FILENAME = "production.xtc"
ATOM_SELECTION = "all"
CUTOFF_NM = 5.0
LAST_FRACTION = 0.25
OUTPUT_PNG = "nearest_neighbors_vs_epsilon.png"


def extract_epsilon(case_dir: Path) -> float:
    match = re.fullmatch(r"pp_eps_(\d+(?:\.\d+)?)", case_dir.name)
    if not match:
        raise ValueError(f"Could not parse epsilon from directory name: {case_dir}")
    return float(match.group(1))


def iter_case_directories(root_dir: Path) -> list[Path]:
    case_dirs = [path for path in root_dir.glob(CASE_GLOB) if path.is_dir()]
    return sorted(case_dirs, key=extract_epsilon)


def start_frame_index(n_frames: int, last_fraction: float) -> int:
    if n_frames <= 0:
        raise ValueError("Trajectory must contain at least one frame")
    if not 0 < last_fraction <= 1:
        raise ValueError("LAST_FRACTION must be greater than 0 and at most 1")
    return int(n_frames * (1.0 - last_fraction))


def cutoff_angstrom(cutoff_nm: float) -> float:
    if cutoff_nm <= 0:
        raise ValueError("CUTOFF_NM must be positive")
    return cutoff_nm * 10.0


def load_dependencies():
    try:
        import matplotlib.pyplot as plt
        import MDAnalysis as mda
        import numpy as np
        from MDAnalysis.lib.distances import self_distance_array
    except ImportError as exc:
        raise SystemExit(
            "This script requires MDAnalysis, matplotlib, and numpy. "
            "Install them before running."
        ) from exc
    return mda, np, plt, self_distance_array


def analyze_case(case_dir: Path, cutoff_nm: float, last_fraction: float, atom_selection: str) -> float:
    mda, np, _, self_distance_array = load_dependencies()

    topology_path = case_dir / TOPOLOGY_FILENAME
    trajectory_path = case_dir / TRAJECTORY_FILENAME
    if not topology_path.exists():
        raise FileNotFoundError(f"Missing topology file: {topology_path}")
    if not trajectory_path.exists():
        raise FileNotFoundError(f"Missing trajectory file: {trajectory_path}")

    universe = mda.Universe(str(topology_path), str(trajectory_path))
    atoms = universe.select_atoms(atom_selection)
    if len(atoms) == 0:
        raise ValueError(f"Selection {atom_selection!r} matched no atoms in {case_dir}")

    trajectory = universe.trajectory
    start = start_frame_index(len(trajectory), last_fraction)
    cutoff_a = cutoff_angstrom(cutoff_nm)
    frame_means = []
    n_atoms = len(atoms)
    upper_rows, upper_cols = np.triu_indices(n_atoms, k=1)

    for ts in trajectory[start:]:
        distances = self_distance_array(atoms.positions, box=ts.dimensions)
        within_cutoff = distances <= cutoff_a
        pair_rows = upper_rows[within_cutoff]
        pair_cols = upper_cols[within_cutoff]
        per_bead_counts = np.bincount(
            np.concatenate((pair_rows, pair_cols)),
            minlength=n_atoms,
        )
        frame_means.append(float(per_bead_counts.mean()))

    if not frame_means:
        raise ValueError(f"No frames selected for analysis in {case_dir}")
    return sum(frame_means) / len(frame_means)


def plot_results(root_dir: Path, epsilons: list[float], averages: list[float], cutoff_nm: float) -> Path:
    _, _, plt, _ = load_dependencies()

    output_path = root_dir / OUTPUT_PNG
    plt.figure(figsize=(8, 5))
    plt.plot(epsilons, averages, marker="o")
    plt.xlabel("Epsilon (kJ/mol)")
    plt.ylabel("Average nearest neighbors")
    plt.title(f"Nearest neighbors within {cutoff_nm} nm")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    return output_path


def main() -> int:
    case_dirs = iter_case_directories(ROOT_DIR)
    if not case_dirs:
        print(f"No directories matching {CASE_GLOB!r} were found in {ROOT_DIR.resolve()}.", file=sys.stderr)
        return 1

    epsilons = []
    averages = []
    for case_dir in case_dirs:
        epsilon = extract_epsilon(case_dir)
        average_neighbors = analyze_case(case_dir, CUTOFF_NM, LAST_FRACTION, ATOM_SELECTION)
        epsilons.append(epsilon)
        averages.append(average_neighbors)
        print(f"{case_dir.name}: epsilon={epsilon:.3f}, average nearest neighbors={average_neighbors:.6f}")

    output_path = plot_results(ROOT_DIR, epsilons, averages, CUTOFF_NM)
    print(f"Saved plot to {output_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
