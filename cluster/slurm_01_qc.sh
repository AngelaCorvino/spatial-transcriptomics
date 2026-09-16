#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --account=dirbas
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

echo "Start: $(date)"

# sbatch copies this script into its spool; use the submission directory.
ROOT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$ROOT_DIR"

JOB_SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/spatial_qc.${SLURM_JOB_ID:-local}.XXXXXX")"
trap 'rm -rf -- "$JOB_SCRATCH"' EXIT
export TMPDIR="$JOB_SCRATCH"
export MPLCONFIGDIR="$TMPDIR/matplotlib"
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export NUMBA_NUM_THREADS="$OMP_NUM_THREADS"

module load miniconda/3
eval "$(conda shell.bash hook)"
conda activate /home/acorvino/.envs/python-env

python -u scripts/01_qc.py --config configs/cluster.yaml "$@"

conda deactivate

echo "Completion: $(date)"
