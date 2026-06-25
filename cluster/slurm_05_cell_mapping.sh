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

eval "$(conda shell.bash hook)"
conda activate spatial-transcriptomics

python -u scripts/05_cell_mapping.py --config configs/cluster.yaml

conda deactivate

echo "Completion: $(date)"
