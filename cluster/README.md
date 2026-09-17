# Minimal Cluster Validation

Use this first to confirm that the repository imports correctly on SCG before
running the full pipeline.

```bash
cd ~
git clone -b main git@github.com:AngelaCorvino/spatial-transcriptomics.git
cd spatial-transcriptomics

mkdir -p /home/acorvino/.envs
mkdir -p /oak/stanford/groups/dirbas/Angela/thymus_9h/raw_data
mkdir -p /labs/dirbas/acorvino/thymus_9h/processed_data

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
Raw data live at `/oak/stanford/groups/dirbas/Angela/thymus_9h/raw_data`;
processed results go to `/labs/dirbas/acorvino/thymus_9h/processed_data`.
Both are hosted on Oak, but `/labs` is SCG-managed storage and the group's
`/oak/stanford/groups` directory belongs to a separate allocation. See the
[Stanford Oak FAQ](https://docs.oak.stanford.edu/faq/). Run `checkquota` on SCG
to check available capacity; this configuration change does not increase quota.

The SLURM smoke test creates per-job scratch space under
`/tmp/$USER/$SLURM_JOB_ID` and uses that location for Matplotlib cache files.
The `logs/` directory is created before submission because SLURM opens output
and error files before the job script starts. Submit all wrappers from the
repository root: they use `SLURM_SUBMIT_DIR` to locate code and configuration
even when SLURM runs a temporary copy of the wrapper. Logs remain in the
repository's `logs/` directory; analysis outputs use the configured results root.

Do not run the full pipeline jobs until `cluster/test_import.sh` completes
successfully.

## Visium HD notebook analysis as a batch job

### Project storage layout

The raw archives and processed results use separate roots:

```text
/oak/stanford/groups/dirbas/Angela/thymus_9h/raw_data/
├── FD1/binned_outputs.tar.gz
├── FD2/binned_outputs.tar.gz
└── ... FD12/

/labs/dirbas/acorvino/thymus_9h/processed_data/
├── qc/visium_hd_008um/
├── preprocess/
├── integration/
├── spatial/
└── cell_mapping/
```

Raw archives stay in their existing `raw_data/FD1`–`raw_data/FD12` folders.
Confirm they are present before submitting a job. Changing this repository's
config does not move cluster files or relocate previously generated results.

```bash
ls -l /oak/stanford/groups/dirbas/Angela/thymus_9h/raw_data/FD{1,2}/binned_outputs.tar.gz
mkdir -p /labs/dirbas/acorvino/thymus_9h/processed_data
```

The scripts and active notebook both read these paths from `configs/cluster.yaml`.
Rerun notebook setup after changing the config. The repository remains at
`/home/acorvino/spatial-transcriptomics`; extracted job inputs remain in `TMPDIR`.

### Copy existing results once

After updating the repository, copy previously generated results to the new
destination **on the cluster**. Wait for jobs writing either results directory
to finish first. Keep the original directory until the copy has been verified.

```bash
mkdir -p /labs/dirbas/acorvino/thymus_9h/processed_data &&
rsync -rltp --chmod=Dg+s --ignore-existing --partial-dir=.rsync-partial --progress \
  /oak/stanford/groups/dirbas/Angela/thymus_9h/processed_data/ \
  /labs/dirbas/acorvino/thymus_9h/processed_data/
```

Use the full paths shown here, including in a new terminal session. A command
such as `"$previous_results/"` becomes `/` when the variable is unset or empty,
which could scan unrelated directories. If output mentions paths outside the
two project directories above, stop with Ctrl+C and inspect the command before
continuing.

The trailing slashes copy the **contents** of `processed_data`, preserving
`qc/visium_hd_008um/FD1/`, etc. Existing destination files are left in place;
the source is retained. Interrupted transfers keep incomplete files in
`.rsync-partial` for resuming, rather than publishing incomplete QC files.
Check file contents with a checksum dry run:

```bash
rsync -rcn --itemize-changes \
  /oak/stanford/groups/dirbas/Angela/thymus_9h/processed_data/ \
  /labs/dirbas/acorvino/thymus_9h/processed_data/
```

This check should exit successfully and print no file differences. If it lists
files, inspect those missing or differing copies before rerunning QC; do not
delete the original results. Extra files at the destination are allowed.

Copy the complete QC directory, including each mouse's `bin_qc.parquet` and
`qc_metadata.json`. With unchanged raw archives and reference lists, the next
QC run can reuse those tables and regenerate plots without reading the matrices.
Do not add `--force` for a plotting-only rerun. If no previous results exist,
skip the copy and run QC normally.

Remove any old `--output-dir` override from submission commands to use the new
configured default. Local downloaded results and previously generated reports
stay where they are; this update changes cluster output destinations.

### Submit QC

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

If a submitted QC job disappears from `squeue`, inspect both
`logs/qc-JOBID.out` and `logs/qc-JOBID.err`; the job may have completed or failed.
An `OD_Jupyt` entry is the separate interactive Jupyter session.
The wrappers suspend Bash's unset-variable check during Module/Conda activation
and deactivation because those hooks may reference an unset interactive variable
such as `PS1`. Strict variable checking resumes before Python runs, and command
failure checking remains enabled throughout setup and analysis.

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
`/labs/dirbas/acorvino/thymus_9h/processed_data/qc/visium_hd_008um/`.
The script creates missing output directories. Override with `--output-dir PATH`
to use that exact QC directory instead (no `visium_hd_008um` suffix is added).

Each mouse folder contains:

- `bin_qc.parquet`: barcodes, counts, detected genes, MT/RP percentages when
  available, and full-resolution pixel coordinates; no expression matrix.
- `qc_metadata.json`: source archive identity, reference genes, and match counts.
- `qc_summary.csv`: QC distributions and percentiles.
- `candidate_summary.csv`: thresholds, retained bins/transcripts, and retained medians.
- `failure_summary.csv`: independent metric failures and their exclusive overlaps
  for every candidate, including the no-additional-filter baseline.
- `qc_distributions.png`: separate distributions of UMI counts, detected genes,
  and mitochondrial percentage before additional filtering. Counts and genes
  use `log10(value + 1)`; mitochondrial percentage stays on its original scale.
- `spatial_qc.png`: the same three metrics mapped over all input bins. Counts
  and genes use `log1p(value)`; mitochondrial percentage is not log-transformed.
- `candidate_spatial_retention.png`: the existing hypothetical retention maps.

Spatial maps draw filled HD bin footprints instead of tiny scatter markers to
avoid a coarse rendering grid. Barcode rows/columns define the bin layout;
an affine fit to image coordinates preserves its orientation. Counts and QC
thresholds are unchanged, with no smoothing or aggregation. Candidate panels
share the color limits of all input bins within each sample. Excluded bins are
gray, and positions absent from the input remain empty.

`comparisons/FD1_FD2/` (or the selected mouse IDs joined by underscores) contains:

- `cross_mouse_qc_summary.csv` and `cross_mouse_candidate_summary.csv`.
- `cross_mouse_failure_summary.csv`: all per-mouse failure tables combined.
- `cross_mouse_retention.png` and `cross_mouse_spatial_*.png` (moderate threshold,
  up to four mice per page).
- `run.json`: requested inputs and run status; `complete` is written only after
  every requested output succeeds. After a failed run, rerun the same command.

### Independent failures and overlaps

`failure_summary.csv` uses the existing candidate rules, with inclusive passing
boundaries. A bin fails UMI or gene QC when its value is below the minimum; it
fails MT QC when its percentage exceeds the maximum. These are diagnostic
classifications, not applied filters. With the current config, the input is
Space Ranger's filtered 8 µm matrix, before additional analysis filtering.

The `group_type` column distinguishes two kinds of rows:

- `exclusive`: `passes_all`, `umi_only`, `genes_only`, `mt_only`,
  `umi_and_genes_only`, `umi_and_mt_only`, `genes_and_mt_only`, and `all_three`.
  Each bin belongs to exactly one of these eight groups for each candidate.
- `marginal`: `fails_umi`, `fails_genes`, and `fails_mt`, each including bins
  that also fail other criteria. These rows overlap and must not be summed.

Each row records bins and UMIs, their percentages of the full input, the input
totals, and the candidate thresholds. Exclusive rows sum to the input totals;
`passes_all` reconciles with `candidate_summary.csv`. UMI percentages are blank
when the input has zero total UMIs. Interpret each candidate separately.

After copying the updated code to SCG, rerun the all-mouse command above without
`--force`. Valid cached per-bin tables are reused to produce the new summaries
and three-metric maps. This update does not require recomputing count metrics.

### Runtime and reuse

Mice are processed sequentially. Only the selected matrix and positions are
extracted into temporary scratch; images are skipped, counts stay sparse, and
unused top-gene QC rankings are disabled. Figures include all bins and are closed
after saving. No filtering threshold is selected automatically.

Completed per-mouse QC tables are reused if the archive path, size, modification
time, matrix selection, and reference lists match. Changing the comparison from
two mice to all 12 reuses completed mice. Plots and summaries are regenerated.
After a plotting-only update, rerun the same command without `--force` to redraw
the figures from the existing per-bin tables when their cache metadata matches.
Use `--force` to recompute metrics, including after code changes that alter QC
semantics or an archive replacement that preserves its size and timestamp.
`--no-plots` writes tables only; existing plots are not updated in that mode.
Avoid simultaneous runs targeting the same output folder.

Archives remain on Oak. Per-mouse scratch is cleaned after processing. The first
run still needs to read/decompress the source archives; moving from a notebook to
a script alone does not make that I/O faster. Cache reuse avoids it on reruns.
