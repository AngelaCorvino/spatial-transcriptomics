#!/bin/bash
#SBATCH --job-name=cellmap
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --account=your_account

set -euo pipefail

echo "Start: $(date)"

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
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

python -u scripts/05_cell_mapping.py --config configs/cluster.yaml

set +u
conda deactivate
set -u

echo "Completion: $(date)"
