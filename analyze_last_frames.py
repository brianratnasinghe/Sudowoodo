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
DEFAULT_PECTIN_ITP = "toppar_custom/sudowoodo_pectin.itp"

# Colour map: MDAnalysis atom type → colour.
# Pectin atomtypes are PctRep, PctNeu, PctXlk (see build_sweep.py).
TYPE_COLOURS = {
    "PctRep": "#e6194b",   # red
    "PctNeu": "#3cb44b",   # green
    "PctXlk": "#4363d8",   # blue
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


def _parse_itp_atomtypes(itp_path: Path) -> dict[str, str]:
    """Parse a GROMACS ITP and return {atom_name: atomtype} from the [atoms] section.

    Falls back to an empty dict if the file cannot be read or has no [atoms] section.
    """
    name_to_type: dict[str, str] = {}
    try:
        in_atoms = False
        for raw_line in itp_path.read_text().splitlines():
            line = raw_line.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("["):
                in_atoms = "atoms" in line.lower()
                continue
            if in_atoms:
                parts = line.split()
                # ITP atoms line: id  type  resnr  residu  atom  cgnr  charge
                if len(parts) >= 5:
                    atomtype = parts[1]
                    atom_name = parts[4]
                    name_to_type[atom_name] = atomtype
    except OSError:
        pass
    return name_to_type


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


def _find_pectin_itp(case_dir: Path) -> Path:
    """Return a path to a pectin ITP in *case_dir*.

    Prefers the monolithic template ``sudowoodo_pectin.itp``; falls back to
    the first per-fiber ITP (``sudowoodo_pectin_1.itp``) produced by the new
    build pipeline.  Returns the template path (possibly non-existent) when
    neither is found so the caller can handle the missing-file case.
    """
    toppar = case_dir / "toppar_custom"
    template = toppar / "sudowoodo_pectin.itp"
    if template.exists():
        return template
    per_fiber = toppar / "sudowoodo_pectin_1.itp"
    if per_fiber.exists():
        return per_fiber
    return template  # let _parse_itp_atomtypes handle the missing file gracefully


def _process_case(case_dir: Path, topo: str, traj: str, selection: str, mda):
    """Load last frame, return (positions, atomtypes, box, case_name)."""
    u = mda.Universe(str(case_dir / topo), str(case_dir / traj))
    # Jump to the last frame
    u.trajectory[-1]
    ag = u.select_atoms(selection)
    positions = ag.positions.copy()   # shape (N, 3) in Angstrom

    # Resolve per-bead atomtypes from the pectin ITP (atom name → atomtype).
    # GRO files don't carry atomtype info, so ag.types would just return the
    # element symbol ("P") for all beads.  The ITP stores the correct types
    # (PctRep / PctNeu / PctXlk).
    itp_path = _find_pectin_itp(case_dir)
    name_to_type = _parse_itp_atomtypes(itp_path)
    if name_to_type:
        atomtypes = [name_to_type.get(name, name) for name in ag.names]
    else:
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
