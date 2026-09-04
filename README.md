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
cp local_config.yaml.template local_config.yaml
# Edit with your local absolute paths
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

- Local execution: `python -u scripts/01_qc.py --config configs/local.yaml`
- Cluster execution: `sbatch cluster/slurm_01_qc.sh`

Local and cluster configs differ in environment settings, output directories,
and cluster-specific submission options, while the analysis code remains shared
in `src/spatial_transcriptomics`.

#### Local vs Cluster config

- `configs/local.yaml` is for running on a laptop or workstation.
- `configs/cluster.yaml` is for running inside an HPC SLURM job.
- Both files define data paths, sample metadata, analysis thresholds, and
  output directories.

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
