# spatial-transcriptomics

A modular Python package for spatial transcriptomics analysis with Jupyter notebook support.

The repository is intentionally kept simple:
- core code in `src/`
- notebooks in `notebooks/`
- tests in `tests/`
- straightforward local commands through `make` (for example:
  `make setup`, `make dev-install`, `make lab`, `make test`, `make lint`,
  `make typecheck`, `make check`)

## Project Structure

```
src/spatial_transcriptomics/       → Reusable analysis modules
├── config.py                       → Path and config management
├── data.py                         → Data loading and preprocessing
├── analysis.py                     → Analysis helpers
└── plotting.py                     → Plotting helpers

notebooks/                          → Analysis workflows
└── visium_thymus_flash_analysis_pipeline.ipynb

tests/                              → Unit tests
```

## Choose Your Setup

Pick one path based on your role. Do not run both.

### End User Setup (run analysis, not developing package code)

1. Clone and enter the repository:

```bash
git clone https://github.com/angelacorvino/spatial-transcriptomics.git
cd spatial-transcriptomics
```

2. Install Python 3.12, then runtime dependencies and the package:

```bash
make setup
```

3. Create your local machine-specific config:

```bash
cp -n configs/local.yaml local_config.yaml
# Edit with your external data/results paths; preserve existing local settings.
```

4. Start Jupyter Lab:

```bash
make lab
```

### Developer Setup (contributing code, tests, lint/type checks)

1. Clone and enter the repository:

```bash
git clone https://github.com/angelacorvino/spatial-transcriptomics.git
cd spatial-transcriptomics
```

2. Install Python 3.12, then runtime + development dependencies:

```bash
make dev-install
```

`make dev-install` includes everything in `make setup` and also installs dev
tools (`pytest`, `ruff`, `mypy`, `jupyterlab`, etc.), so developers should run
`make dev-install` only.

Local Make commands use `.venv-dev/` and Python 3.12, matching CircleCI.
The existing `.venv/` is retained; use `make VENV=.venv lab` to open it.
Development and CI installs share `constraints-dev.txt`, including transitive
dependencies. These constraints do not apply to the cluster's Conda environment
or its `pip install -e .` workflow. Python patch versions and platform-specific
packages may differ between macOS and Linux.

To update dependency pins, install `pip setuptools wheel '.[dev]'` in a fresh
Python 3.12 environment. Run `make check` with `PYTHON` pointing to that Python, then
regenerate the pins with `python -m pip freeze --all --exclude-editable` (retain
the constraints file's header). Verify a constrained install and the CI checks
before adopting the update.

3. Create your local machine-specific config:

```bash
cp local_config.yaml.template local_config.yaml
# Edit with your local absolute paths
```

4. Run checks during development:

```bash
make test
make lint
make typecheck
# or all at once:
make check
```

5. Start Jupyter Lab for notebook work:

```bash
make lab
```

## Workflow

### Notebook Workflow

Notebooks are for exploration, QC inspection, and plotting only. They should
call reusable functions from `src/spatial_transcriptomics` and not duplicate
analysis logic found in the `scripts/` pipeline.

```python
from spatial_transcriptomics.config import load_config
from spatial_transcriptomics.data import load_data, preprocess_data
from spatial_transcriptomics.analysis import compute_statistics

config = load_config()
```

Complete setup first using the appropriate path in **Choose Your Setup**, then
start Jupyter from the repository root:

```bash
make lab
```

### Script-based Reproducible Workflow

This repository supports both local and cluster execution without duplicating
analysis code.

- Local execution: `python -u scripts/01_qc.py --config local_config.yaml`
- Cluster Visium HD QC: `sbatch cluster/slurm_01_qc.sh --visium-hd --mice FD1 FD2`

To review moderate QC failures against H&E, add `--he-overlays`.
Each sample gets `qc_failures_he.png`: H&E alone, MT-only, low-content only
(UMI < 25 or genes < 20), and overlapping failures (MT > 20%). No bins are removed.
The script reuses cached QC and reads the 8-µm scale factor from
`binned_outputs.tar.gz`. If `spatial.tar.gz` is stored elsewhere, specify its
parent sample directory root:

```bash
sbatch cluster/slurm_01_qc.sh --visium-hd --mice FD5 FD9 --he-overlays \
  --he-image-root /oak/stanford/groups/dirbas/raw/Dirbas_Processed_Visium_HD_Data
```

Check tissue edges and internal landmarks for alignment before interpreting
the overlays. These figures show the filter's expected effects without applying it.

#### Apply the selected thymus QC filter

QC filter v1 (selected 19 October 2026) retains Space Ranger filtered, in-tissue
8-µm bins with **UMIs ≥25, detected genes ≥20, and MT ≤20%**. All three boundaries
are inclusive. It keeps all tissue eligible, including pink tissue, edges and
uncertain fragments; proposed anatomical contours remain review-only. There is
no upper UMI/gene cutoff, ribosomal cutoff, or whole-sample exclusion. A QC failure
is an operational exclusion, not proof of tissue damage.

Run on the cluster after the existing HD QC tables are available:

```bash
sbatch cluster/slurm_02_preprocess.sh --visium-hd --all-mice
```

This HD branch uses the selected v1 rules rather than the generic `analysis`
defaults. It reloads one mouse at a time, verifies barcodes, metrics and spatial
coordinates against the saved QC review, and preserves the inputs. Float64 MT
arithmetic retains exact-20% bins that older rounding could exclude; the summary
records these boundary corrections. Other changes in retention stop the run.
Outputs are under the configured `paths.preprocess_output`, in
`visium_hd_008um/thymus_8um_qc_v1/`:

- `FD*/filtered_counts.h5ad`: retained **unnormalized** counts, all gene features,
  sample IDs, full-resolution coordinates and image scale factors.
- `FD*/bin_selection.parquet`: metrics and independent failure flags for every
  input bin, including excluded bins.
- `filter_summary.csv` and `FD*/failure_summary.csv`: per-mouse effects.
- `run.json`: filter definition, input paths, code hashes, execution status and
  completed mice. A complete 12-mouse run is
  required before treating application as finished.

Existing output directories are preserved. For a retry, supply an unused
`--output-dir`; `--source-root` and `--qc-dir` can locate relocated inputs.

### Tissue annotation after Visium HD filtering

The first pass identifies broad cell types from cluster markers and H&E. Published
[Visium HD work](https://www.nature.com/articles/s41588-025-02193-3) used 8-µm
bins and marker-based annotation to identify ten broad cell types.

The existing clustering script accepts the completed QC v1 matrices. It processes
each mouse independently: normalize to 10,000 counts per bin, log2-transform,
select 2,000 variable genes, and run PCA, neighbours, Leiden (resolution 0.5),
and UMAP. MT genes are excluded from the embedding; all retained bins, gene
features, and original counts are preserved. No batch correction is applied.

Install the clustering extra in the cluster environment, then submit from the
repository root:

```bash
conda activate /home/acorvino/.envs/python-env
python -m pip install -e '.[annotation]'
mkdir -p logs
sbatch cluster/slurm_03_integration.sh --visium-hd --all-mice \
  --he-source-root /oak/stanford/groups/dirbas/raw/Dirbas_Processed_Visium_HD_Data
```

Outputs use `integration/visium_hd_008um/FD1` through `FD12`:

- `clustered_bins.h5ad`: log-normalized expression, raw counts in `layers['counts']`,
  spatial coordinates, PCA, graph, UMAP and sample-specific cluster IDs.
- `clusters.png`, `marker_dotplot.png`, `marker_spatial.png`: annotation views.
- `qc_failures_he.png`: actual QC v1 failures on H&E, when `--he-source-root` is given.
- `marker_coverage.csv`, `marker_expression.csv`, `cluster_markers.csv`: marker evidence.
- `cluster_annotation.csv`: cluster sizes, QC metrics and top genes, with broad
  cell-type/tissue labels to complete after reviewing the plots and markers.

The root contains `annotation_summary.csv` and `run.json`. Download the PNG/CSV
files for review; matrices can remain on the cluster. The workflow refuses to
overwrite an existing output directory; use `--output-dir` for a separate rerun.

[Marker panels](configs/thymus_markers.yaml) combine the supplied lung marker tables
with mouse thymus references. Version `thymus_annotation_markers_v1` records each
gene's source, source-table hashes and interpretation notes; the complete selection
is captured in `run.json`. The panel guides annotation plots, not clustering.
Missing genes are recorded separately from present genes with no detected counts.
Assign broad labels using multiple markers, spatial location and H&E; retain a
mixed or uncertain label where needed. Cluster IDs differ between mice.
Marker ranks describe cluster-vs-rest expression, not treatment effects.

Implementation follows [Scanpy's clustering workflow](https://scanpy.scverse.org/en/stable/tutorials/basics/clustering.html).
Use `--mice FD5 FD9` for a subset; resolution and embedding dimensions can be
changed through the CLI without changing the QC filter.

Local and cluster configs differ in environment settings, output directories,
and cluster-specific submission options, while the analysis code remains shared
in `src/spatial_transcriptomics`.

#### Local vs Cluster config

- `configs/local.yaml` is a laptop example. For this project, copy it to the
  gitignored `local_config.yaml` and set `data_dir`, `output_dir`, `visium_hd`
  staging/source paths, and `paths` outputs outside the code repository.
- Set `report_dir` in `local_config.yaml` to the experimental report folder so
  the report agent can locate it. Keep reports and their images with the results.
- `configs/cluster.yaml` is for running inside an HPC SLURM job.
- Both files define data paths, sample metadata, analysis thresholds, and
  output directories.

#### Thymus project storage on SCG

`configs/cluster.yaml` supplies the shared paths for the scripts and the active
Visium HD notebook:

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

Both locations are hosted on Oak: `/labs` is SCG-managed storage, while the
`/oak/stanford/groups` path belongs to a separate Oak allocation. See the
[Stanford Oak FAQ](https://docs.oak.stanford.edu/faq/). Check available quota on
the cluster; changing paths does not establish that more space is available.

FD1–FD12 sample folders remain under the raw-data path above. QC reports are
saved under the `/labs` results root in `qc/visium_hd_008um`; the other output
directories are reserved for their corresponding pipeline steps.
Updating the config does not move existing results: copy them as described in
the cluster migration instructions to reuse cached QC tables.
Extraction uses job-local `TMPDIR`. Repository code and the Python environment
remain under `/home/acorvino/`.

See [cluster/README.md](cluster/README.md) for folder setup, migration notes, and
job submission commands. Create `logs/` in the repository before submitting jobs.

Because notebooks are exploratory, reproducible work should use the scripts in
`scripts/` and the configs in `configs/`.

Because the package is installed in editable mode, notebook code can import
`spatial_transcriptomics` directly without modifying `sys.path`.

See [notebooks/README.md](notebooks/README.md) for detailed workflow documentation.

### Commit Workflow

From the repository root:

```bash
make check
git status
git add .
git commit -m "Describe the analysis change"
git push
```

If you changed notebook outputs only, clear them before committing to keep diffs
small.

## Testing

```bash
make test
```

## Linting and Type Checking

```bash
make lint
make typecheck
```

Both commands target `src/` and `tests/`, which matches the CI pipeline.
