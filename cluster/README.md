# Minimal Cluster Validation

Use this first to confirm that the repository imports correctly on SCG before
running the full pipeline.

```bash
cd ~
git clone git@github.com:YOUR_USERNAME/spatial-transcriptomics.git
cd spatial-transcriptomics

mkdir -p /labs/dirbas/$USER/.envs
mkdir -p /labs/dirbas/$USER/spatial_data
mkdir -p /labs/dirbas/$USER/spatial_results

conda env create -p /labs/dirbas/$USER/.envs/spatial-transcriptomics -f environment.yml
conda activate /labs/dirbas/$USER/.envs/spatial-transcriptomics
pip install -e .

python scripts/test_import.py
mkdir -p logs
sbatch cluster/test_import.sh
squeue -u $USER
ls logs
```

Repository code should live under `/home/$USER/spatial-transcriptomics`.
Conda environments, large data, and results should live under
`/labs/dirbas/$USER/`.

The SLURM smoke test creates per-job scratch space under
`/tmp/$USER/$SLURM_JOB_ID` and uses that location for Matplotlib cache files.
The `logs/` directory is created before submission because SLURM opens output
and error files before the job script starts.

Do not run the full pipeline jobs until `cluster/test_import.sh` completes
successfully.
