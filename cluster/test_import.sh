#!/bin/bash
#SBATCH --job-name=st_test
#SBATCH --partition=batch
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --account=dirbas
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

mkdir -p logs

export TMPDIR="/tmp/$USER/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

export MPLCONFIGDIR="$TMPDIR/matplotlib"
mkdir -p "$MPLCONFIGDIR"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export NUMBA_NUM_THREADS="$SLURM_CPUS_PER_TASK"

echo "Job started at $(date)"
echo "Repository root: $ROOT_DIR"
echo "Temporary directory: $TMPDIR"

module load miniconda/3
eval "$(conda shell.bash hook)"
conda activate "/labs/dirbas/$USER/.envs/spatial-transcriptomics"
echo "Conda environment: $CONDA_PREFIX"

python -u scripts/test_import.py

conda deactivate

echo "Job finished at $(date)"
