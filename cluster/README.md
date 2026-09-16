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

## Visium HD notebook analysis as a batch job

The existing `scripts/01_qc.py --visium-hd` reproduces the 8 µm notebook QC
and candidate-threshold comparison without a Jupyter session. It does not apply
filters, normalize counts, or cluster bins. The generic QC script behavior remains
available when `--visium-hd` is omitted.

From the repository root **on SCG**, after copying the updated repository code:

```bash
cd /home/acorvino/spatial-transcriptomics
mkdir -p logs
sbatch cluster/slurm_01_qc.sh --visium-hd --mice FD1 FD2
```

Once the two-mouse comparison has been reviewed, run all 12:

```bash
sbatch cluster/slurm_01_qc.sh --visium-hd --all-mice
```

The job requests 4 CPUs, 32 GB RAM, and 2 hours using account `dirbas`.
These are starting allocations, not measured requirements for all 12 samples;
override with `sbatch --mem=64G --time=04:00:00 ...` if needed. Submit from the
repository root because the wrapper uses `SLURM_SUBMIT_DIR`. Follow the job with
`squeue -u "$USER"` and `tail -f logs/qc-JOBID.out` (substitute its job ID).

Alternatively, run in the terminal of an **existing Jupyter compute session**:

```bash
cd /home/acorvino/spatial-transcriptomics
/home/acorvino/.envs/python-env/bin/python -u scripts/01_qc.py \
  --config configs/cluster.yaml --visium-hd --mice FD1 FD2
```

That command uses the session's `TMPDIR` and must finish before its allocation
ends. Use the batch submission to run independently of Jupyter. Do not run the
analysis on a login node.

### Saved outputs

By default, outputs go to
`/oak/stanford/groups/dirbas/acorvino/spatial_results/qc/visium_hd_008um/`.
Override with `--output-dir PATH`.

Each mouse folder contains:

- `bin_qc.parquet`: barcodes, counts, detected genes, MT/RP percentages when
  available, and full-resolution pixel coordinates; no expression matrix.
- `qc_metadata.json`: source archive identity, reference genes, and match counts.
- `qc_summary.csv`: QC distributions and percentiles.
- `candidate_summary.csv`: thresholds, retained bins/transcripts, and retained medians.
- `qc_distributions.png`, `spatial_qc.png`, and `candidate_spatial_retention.png`.

`comparisons/FD1_FD2/` (or the selected mouse IDs joined by underscores) contains:

- `cross_mouse_qc_summary.csv` and `cross_mouse_candidate_summary.csv`.
- `cross_mouse_retention.png` and `cross_mouse_spatial_*.png` (moderate threshold,
  up to four mice per page).
- `run.json`: requested inputs and run status; `complete` is written only after
  every requested output succeeds. After a failed run, rerun the same command.

### Runtime and reuse

Mice are processed sequentially. Only the selected matrix and positions are
extracted into temporary scratch; images are skipped, counts stay sparse, and
unused top-gene QC rankings are disabled. Figures include all bins and are closed
after saving. No filtering threshold is selected automatically.

Completed per-mouse QC tables are reused if the archive path, size, modification
time, matrix selection, and reference lists match. Changing the comparison from
two mice to all 12 reuses completed mice. Plots and summaries are regenerated.
Use `--force` to recompute metrics, including after code changes that alter QC
semantics or an archive replacement that preserves its size and timestamp.
`--no-plots` writes tables only; existing plots are not updated in that mode.
Avoid simultaneous runs targeting the same output folder.

Archives remain on Oak. Per-mouse scratch is cleaned after processing. The first
run still needs to read/decompress the source archives; moving from a notebook to
a script alone does not make that I/O faster. Cache reuse avoids it on reruns.
