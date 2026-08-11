#!/bin/bash
# submit_sweep.sh — write a SLURM submit.sh in each sweep subdirectory and submit it.
#
# Usage:
#   bash submit_sweep.sh <sweep_dir> [options]
#
# Arguments:
#   sweep_dir        Parent directory produced by build_bead_count_sweep.py (required)
#
# Options:
#   --partition P    SLURM partition (default: astro2_long)
#   --time T         Wall-clock limit (default: 72:00:00)
#   --cpus N         CPUs per task / OMP threads (default: 24)
#   --dry-run        Write submit.sh files but do not call sbatch
#
# Example:
#   bash submit_sweep.sh my_sweep --partition astro2_long --time 48:00:00 --cpus 24
#   bash submit_sweep.sh my_sweep --dry-run

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
PARTITION="astro2_long"
TIME="72:00:00"
CPUS=24
DRY_RUN=0
SWEEP_DIR=""

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --partition) PARTITION="$2"; shift 2 ;;
        --time)      TIME="$2";      shift 2 ;;
        --cpus)      CPUS="$2";      shift 2 ;;
        --dry-run)   DRY_RUN=1;      shift   ;;
        -*)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
        *)
            if [[ -z "$SWEEP_DIR" ]]; then
                SWEEP_DIR="$1"
                shift
            else
                echo "Unexpected argument: $1" >&2
                exit 1
            fi
            ;;
    esac
done

if [[ -z "$SWEEP_DIR" ]]; then
    echo "Usage: bash submit_sweep.sh <sweep_dir> [--partition P] [--time T] [--cpus N] [--dry-run]" >&2
    exit 1
fi

if [[ ! -d "$SWEEP_DIR" ]]; then
    echo "Error: directory not found: $SWEEP_DIR" >&2
    exit 1
fi

SWEEP_DIR="$(realpath "$SWEEP_DIR")"

# ---------------------------------------------------------------------------
# Submit one job per subdirectory
# ---------------------------------------------------------------------------
submitted=0
skipped=0

for case_dir in "$SWEEP_DIR"/*/; do
    [[ -d "$case_dir" ]] || continue

    run_sh="$case_dir/run.sh"
    if [[ ! -f "$run_sh" ]]; then
        echo "[skip] $case_dir — no run.sh found"
        skipped=$((skipped + 1))
        continue
    fi

    job_name="$(basename "$case_dir")"
    abs_case_dir="$(realpath "$case_dir")"
    submit_sh="$case_dir/submit.sh"

    cat > "$submit_sh" <<SLURM
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --time=${TIME}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --threads-per-core=1
#SBATCH --chdir=${abs_case_dir}
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module purge
module load astro

export SRUN_CPUS_PER_TASK=\$SLURM_CPUS_PER_TASK

source /software/astro/gromacs/2023.3-cpu-no-mpi/bin/GMXRC

bash run.sh
SLURM

    chmod 755 "$submit_sh"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[dry-run] would submit: $submit_sh"
    else
        sbatch "$submit_sh"
        echo "[submitted] $job_name"
        submitted=$((submitted + 1))
    fi
done

echo ""
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[done] dry-run complete — submit.sh written in each case directory"
else
    echo "[done] submitted $submitted jobs, skipped $skipped"
fi
