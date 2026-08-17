#!/usr/bin/env python3
"""analyze_last_frames.py — render the last production frame (pectin only) for
every sweep case directory and save each frame as an individual PNG file.

Each PNG shows a 2-D (X-Y) projection of the pectin bead positions in the
last frame of that simulation, coloured by bead type (PctRep / PctNeu / PctXlk).

Usage
-----
    python analyze_last_frames.py --sweep <sweep_dir> [options]

Arguments
---------
    --sweep DIR     Parent directory produced by build_bead_count_sweep.py or
                    afm_build_sweep.py.  The script searches one level deep for
                    subdirectories that contain afm_system.gro + production.xtc.
    --out-dir DIR   Output directory for PNG files (default: last_frames/)
    --topo FILE     Topology filename to look for (default: afm_system.gro)
    --traj FILE     Trajectory filename to look for (default: production.xtc)
    --selection SEL MDAnalysis selection for pectin beads
                    (default: "resname Pctn")
    --marker-size N Scatter plot marker size (default: 6)
    --dpi N         PNG resolution in dots per inch (default: 150)

Requirements
------------
    MDAnalysis, matplotlib  (pip install MDAnalysis matplotlib)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_TOPO = "afm_system.gro"
DEFAULT_TRAJ = "production.xtc"
DEFAULT_SELECTION = "resname Pctn"
DEFAULT_OUT_DIR = "last_frames"
DEFAULT_MARKER_SIZE = 6
DEFAULT_DPI = 150

# Colour map: MDAnalysis atom type → colour.
# Pectin atomtypes are PctRep, PctNeu, PctXlk (see build_sweep.py).
TYPE_COLOURS = {
    "PctRep": "#e6194b",   # red
    "PctNeu": "#4363d8",   # blue
    "PctXlk": "#3cb44b",   # green
}
DEFAULT_COLOUR = "#aaaaaa"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_deps():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit(
            "This script requires MDAnalysis and matplotlib.\n"
            "Install with:  pip install MDAnalysis matplotlib"
        ) from exc
    return mda, plt


def _find_cases(sweep_dir: Path, topo: str, traj: str) -> list[Path]:
    """Return sorted list of subdirectories that contain topo + traj files."""
    candidates = []
    for sub in sorted(sweep_dir.iterdir()):
        if not sub.is_dir():
            continue
        if (sub / topo).exists() and (sub / traj).exists():
            candidates.append(sub)
    return candidates


def _atom_colour(atomtype: str) -> str:
    return TYPE_COLOURS.get(atomtype, DEFAULT_COLOUR)


def _render_frame(ax, positions, atomtypes, box, case_name: str) -> None:
    """Draw a 2-D X-Y scatter of pectin bead positions onto *ax*."""
    if len(positions) == 0:
        ax.text(0.5, 0.5, "no pectin atoms found",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title(case_name)
        return

    colours = [_atom_colour(t) for t in atomtypes]
    ax.scatter(positions[:, 0], positions[:, 1],
               c=colours, s=DEFAULT_MARKER_SIZE, linewidths=0)

    # Draw box outline
    ax.set_xlim(0, box[0])
    ax.set_ylim(0, box[1])

    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_title(case_name, fontsize=9)
    ax.set_aspect("equal")

    # Legend
    seen_types = dict.fromkeys(atomtypes)  # preserves insertion order
    handles = [
        __import__("matplotlib.patches", fromlist=["Patch"]).Patch(
            facecolor=_atom_colour(t), label=t
        )
        for t in seen_types
    ]
    if handles:
        ax.legend(handles=handles, fontsize=7, loc="upper right",
                  markerscale=1, framealpha=0.7)


def _process_case(case_dir: Path, topo: str, traj: str, selection: str, mda):
    """Load last frame, return (positions, atomtypes, box, case_name)."""
    u = mda.Universe(str(case_dir / topo), str(case_dir / traj))
    # Jump to the last frame
    u.trajectory[-1]
    ag = u.select_atoms(selection)
    positions = ag.positions.copy()   # shape (N, 3) in Angstrom
    atomtypes = list(ag.types)
    box = u.dimensions[:3]            # (Lx, Ly, Lz) in Angstrom
    return positions, atomtypes, box


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Render the last production frame (pectin only) for every sweep "
            "case directory and save each frame as an individual PNG file."
        )
    )
    p.add_argument("--sweep", type=Path, required=True,
                   help="Parent sweep directory to search for cases")
    p.add_argument("--out-dir", type=Path, default=Path(DEFAULT_OUT_DIR),
                   help=f"Output directory for PNG files (default: {DEFAULT_OUT_DIR}/)")
    p.add_argument("--topo", default=DEFAULT_TOPO,
                   help=f"Topology filename (default: {DEFAULT_TOPO})")
    p.add_argument("--traj", default=DEFAULT_TRAJ,
                   help=f"Trajectory filename (default: {DEFAULT_TRAJ})")
    p.add_argument("--selection", default=DEFAULT_SELECTION,
                   help=f"MDAnalysis atom selection (default: '{DEFAULT_SELECTION}')")
    p.add_argument("--marker-size", type=float, default=DEFAULT_MARKER_SIZE,
                   help=f"Scatter marker size (default: {DEFAULT_MARKER_SIZE})")
    p.add_argument("--dpi", type=int, default=DEFAULT_DPI,
                   help=f"PNG resolution in dots per inch (default: {DEFAULT_DPI})")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.sweep.is_dir():
        sys.exit(f"Error: sweep directory not found: {args.sweep}")

    mda, plt = _load_deps()

    cases = _find_cases(args.sweep, args.topo, args.traj)
    if not cases:
        sys.exit(
            f"No case directories found under {args.sweep} "
            f"with both {args.topo} and {args.traj}."
        )

    print(f"[INFO] Found {len(cases)} case(s) in {args.sweep}")

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    global DEFAULT_MARKER_SIZE
    DEFAULT_MARKER_SIZE = args.marker_size

    saved = 0
    for case_dir in cases:
        case_name = case_dir.name
        print(f"  Processing {case_name} ...", end=" ", flush=True)
        try:
            positions, atomtypes, box = _process_case(
                case_dir, args.topo, args.traj, args.selection, mda
            )
        except Exception as exc:
            print(f"SKIPPED ({exc})")
            continue

        fig, ax = plt.subplots(figsize=(8, 6))
        _render_frame(ax, positions, atomtypes, box, case_name)
        fig.tight_layout()
        out_path = out_dir / f"{case_name}.png"
        fig.savefig(str(out_path), dpi=args.dpi)
        plt.close(fig)
        saved += 1
        print(f"ok  ({len(positions)} pectin beads) → {out_path}")

    print(f"\n[DONE] {saved} PNG(s) written to {out_dir}/")


if __name__ == "__main__":
    main()
