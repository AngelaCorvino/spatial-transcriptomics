#!/bin/bash
#SBATCH --job-name=integration
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --account=your_account

set -euo pipefail

echo "Start: $(date)"

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

module load miniconda/3
eval "$(conda shell.bash hook)"
conda activate /home/acorvino/.envs/python-env

python -u scripts/03_integration.py --config configs/cluster.yaml

conda deactivate

echo "Completion: $(date)"
