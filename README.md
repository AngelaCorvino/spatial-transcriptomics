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

2. Install runtime dependencies and the package:

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

2. Install runtime + development dependencies:

```bash
make dev-install
```

`make dev-install` includes everything in `make setup` and also installs dev
tools (`pytest`, `ruff`, `mypy`, `jupyterlab`, etc.), so developers should run
`make dev-install` only.

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
