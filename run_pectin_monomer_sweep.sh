#!/bin/bash
# run_pectin_monomer_sweep.sh
# Sequentially run the GROMACS simulation in each pectin monomer sweep case directory.
#
# Usage: bash run_pectin_monomer_sweep.sh [sweep_dir]
#   sweep_dir  Parent directory containing the case subdirectories
#              (default: pectin_monomer_sweep)

set -euo pipefail

SWEEP_DIR="${1:-pectin_monomer_sweep}"

if [[ ! -d "$SWEEP_DIR" ]]; then
    echo "[error] Sweep directory not found: $SWEEP_DIR" >&2
    exit 1
fi

# Collect case subdirectories that contain a run.sh, in sorted order
mapfile -t CASES < <(
    find "$SWEEP_DIR" -mindepth 1 -maxdepth 1 -type d | sort
)

if [[ ${#CASES[@]} -eq 0 ]]; then
    echo "[error] No subdirectories found in $SWEEP_DIR" >&2
    exit 1
fi

echo "[info] Found ${#CASES[@]} case(s) in $SWEEP_DIR"

HAS_FAILED=0
for CASE_DIR in "${CASES[@]}"; do
    if [[ ! -f "$CASE_DIR/run.sh" ]]; then
        echo "[skip] No run.sh in $CASE_DIR — skipping"
        continue
    fi

    echo "[run ] $(date '+%Y-%m-%d %H:%M:%S')  $CASE_DIR"
    if (cd "$CASE_DIR" && bash run.sh); then
        echo "[done] $(date '+%Y-%m-%d %H:%M:%S')  $CASE_DIR"
    else
        echo "[fail] $(date '+%Y-%m-%d %H:%M:%S')  $CASE_DIR — non-zero exit, stopping." >&2
        HAS_FAILED=1
        break
    fi
done

if [[ $HAS_FAILED -eq 0 ]]; then
    echo "[info] All cases completed successfully."
else
    exit 1
fi
