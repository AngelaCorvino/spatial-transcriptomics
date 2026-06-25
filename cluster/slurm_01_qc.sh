#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=your_account

set -euo pipefail

echo "Start: $(date)"

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

eval "$(conda shell.bash hook)"
conda activate spatial-transcriptomics

python -u scripts/01_qc.py --config configs/cluster.yaml

conda deactivate

echo "Completion: $(date)"
