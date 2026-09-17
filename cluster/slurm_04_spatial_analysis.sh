#!/bin/bash
#SBATCH --job-name=spatial
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

# Module/Conda initialization may read unset interactive variables such as PS1.
# Keep error checking enabled, but suspend nounset while running these hooks.
set +u
module load miniconda/3
CONDA_BASH_HOOK="$(conda shell.bash hook)"
eval "$CONDA_BASH_HOOK"
unset CONDA_BASH_HOOK
conda activate /home/acorvino/.envs/python-env
set -u

python -u scripts/04_spatial_analysis.py --config configs/cluster.yaml

set +u
conda deactivate
set -u

echo "Completion: $(date)"
