# spatial-transcriptomics

A modular Python package for spatial transcriptomics analysis with Jupyter notebook support.

## Project Structure

```
src/spatial_transcriptomics/       → Reusable analysis modules
├── config.py                       → Path and config management
├── data.py                         → Data loading and preprocessing
├── analysis.py                     → Analysis utilities
└── example.py                      → Example utility functions

notebooks/                          → Analysis workflows
└── visium_thymus_flash_analysis_pipeline.ipynb

tests/                              → Unit tests
```

## Installation

```bash
make setup
```

## Development

```bash
make dev-install
```

## Workflow

### Simple Notebook Workflow

```python
from spatial_transcriptomics.config import load_config
from spatial_transcriptomics.data import load_data, preprocess_data
from spatial_transcriptomics.analysis import compute_statistics

config = load_config()
```

Use this once per machine:

```bash
make dev-install
```

This creates a local `.venv/`, upgrades the Python packaging tools, and installs
the package in editable mode.

Then start Jupyter from the repository root:

```bash
make lab
```

Because the package is installed in editable mode, notebook code can import
`spatial_transcriptomics` directly without modifying `sys.path`.

### Local Configuration

Create `local_config.yaml` for absolute paths (gitignored):

```bash
cp local_config.yaml.template local_config.yaml
# Edit with your local absolute paths
```

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
