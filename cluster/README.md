# Minimal Cluster Validation

Use this first to confirm that the repository imports correctly on SCG before
running the full pipeline.

```bash
cd ~
git clone -b cluster-test git@github.com:AngelaCorvino/spatial-transcriptomics.git
cd spatial-transcriptomics

mkdir -p /home/acorvino/.envs
mkdir -p /oak/stanford/groups/dirbas/acorvino/spatial_data
mkdir -p /oak/stanford/groups/dirbas/acorvino/spatial_results

module load miniconda/3
eval "$(conda shell.bash hook)"
conda create -p /home/acorvino/.envs/python-env python=3.11 pip -y
conda activate /home/acorvino/.envs/python-env

python -m pip install --upgrade pip setuptools wheel
python -m pip install -e .
python scripts/test_import.py
mkdir -p logs
sbatch cluster/test_import.sh
squeue -u $USER
ls logs
```

Repository code should live under `/home/acorvino/spatial-transcriptomics`.
The conda environment should live under `/home/acorvino/.envs/python-env`.
Large data and results should live on Oak:
`/oak/stanford/groups/dirbas/acorvino/spatial_data` and
`/oak/stanford/groups/dirbas/acorvino/spatial_results`.

The SLURM smoke test creates per-job scratch space under
`/tmp/$USER/$SLURM_JOB_ID` and uses that location for Matplotlib cache files.
The `logs/` directory is created before submission because SLURM opens output
and error files before the job script starts.

Do not run the full pipeline jobs until `cluster/test_import.sh` completes
successfully.
