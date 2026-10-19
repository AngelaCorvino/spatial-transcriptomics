#!/bin/bash
#SBATCH --job-name=integration
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --account=dirbas
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

echo "Start: $(date)"

# Submit from the repository root; SLURM runs a spooled copy of this script.
ROOT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$ROOT_DIR"

JOB_SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/spatial_annotation.${SLURM_JOB_ID:-local}.XXXXXX")"
trap 'rm -rf -- "$JOB_SCRATCH"' EXIT
export TMPDIR="$JOB_SCRATCH"
export MPLCONFIGDIR="$TMPDIR/matplotlib"
export NUMBA_CACHE_DIR="$TMPDIR/numba"
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export NUMBA_NUM_THREADS="$OMP_NUM_THREADS"

# Module/Conda initialization may read unset interactive variables such as PS1.
# Keep error checking enabled, but suspend nounset while running these hooks.
set +u
module load miniconda/3
CONDA_BASH_HOOK="$(conda shell.bash hook)"
eval "$CONDA_BASH_HOOK"
unset CONDA_BASH_HOOK
conda activate /home/acorvino/.envs/python-env
set -u

python -u scripts/03_integration.py --config configs/cluster.yaml "$@"

set +u
conda deactivate
set -u

echo "Completion: $(date)"
