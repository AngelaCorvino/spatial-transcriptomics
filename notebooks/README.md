# Notebooks

This folder contains Jupyter notebooks for analysis and exploration using the `spatial_transcriptomics` package.

## Workflow

All analysis notebooks should import functions from the installed
`spatial_transcriptomics` package instead of defining analysis logic inline.

### Typical Usage

```python
from spatial_transcriptomics.config import load_config
from spatial_transcriptomics.data import load_data, preprocess_data
from spatial_transcriptomics.analysis import compute_statistics

config = load_config()
data = load_data(config['data_dir'] + '/sample.h5ad')
adata = preprocess_data(data)
stats = compute_statistics(adata)
```

## Configuration

- **No hardcoded paths**: All paths are managed via `local_config.yaml` or defaults
- **Local paths only**: Copy `../local_config.yaml.template` to `../local_config.yaml` and add your absolute paths
- **Gitignored**: `local_config.yaml` is gitignored—safe for local-only absolute paths

### Setup

```bash
cd ..
make dev-install
cp ../local_config.yaml.template ../local_config.yaml
make lab
```

Then open the notebook from Jupyter Lab. Imports should work directly because
the package is installed in editable mode.

## Notebook Best Practices

1. **Import from package**: Keep reusable logic in `src/spatial_transcriptomics/`
2. **Use config for paths**: Never hardcode absolute paths; use `load_config()` or `get_data_path()`
3. **Clear outputs before commit**: Remove execution results to keep notebook files small
4. **Run checks before commit**: From the repo root, run `make check`
5. **Document assumptions**: Add markdown cells explaining data format and dependencies

## Available Functions

| Module | Function | Purpose |
|--------|----------|---------|
| `config` | `load_config()` | Load paths and settings |
| `config` | `get_data_path(filename)` | Get full path to data file |
| `data` | `load_data(path)` | Load spatial data |
| `data` | `preprocess_data(data)` | Run preprocessing pipeline |
| `data` | `filter_by_quality(data, min_counts)` | QC filtering |
| `analysis` | `compute_statistics(data)` | Compute summary stats |
| `analysis` | `prepare_results(data, stats)` | Format results for export |

See `src/spatial_transcriptomics/` for full documentation.
