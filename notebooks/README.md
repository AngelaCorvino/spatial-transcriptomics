# Notebooks

This folder contains Jupyter notebooks for analysis and exploration using the `spatial_transcriptomics` package.

## Workflow

All analysis notebooks should import functions from the installed
`spatial_transcriptomics` package instead of defining analysis logic inline.

### Generic local example

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

- **Cluster notebook**: `visium_thymus_flash_analysis_pipeline.ipynb` defaults to
  `configs/cluster.yaml`. Set `SPATIAL_CONFIG` before starting Jupyter to override it.
- **Raw data**: `/oak/stanford/groups/dirbas/Angela/thymus_9h/raw_data/`, with each
  mouse archive at `FD1/binned_outputs.tar.gz`, `FD2/binned_outputs.tar.gz`, etc.
- **Processed data**: `/labs/dirbas/acorvino/thymus_9h/processed_data/`.
  The batch QC script saves reports under `qc/visium_hd_008um/`; the notebook
  displays its QC tables and plots without automatically exporting them.
- **Scratch**: `${TMPDIR}/FD1_HD` holds the staged FD1 inputs for the current job.
- **Local use**: `configs/local.yaml` provides laptop settings. Calling
  `load_config()` without an argument instead reads the gitignored
  `local_config.yaml`, or repository defaults when that file is absent.

Raw archives remain at the path above. Both `/labs` and the group's `/oak`
allocation are hosted on Oak. Rerun the notebook setup cell after updating the
configuration so the kernel uses the new paths. Existing batch results must be
copied separately to the new output directory to reuse QC caches; see
[cluster setup](../cluster/README.md) for migration and batch jobs.

### Local setup (from the repository root)

```bash
make dev-install
# First-time setup only; preserve an existing local_config.yaml.
cp -n configs/local.yaml local_config.yaml
# Set data, output, and staging paths outside the code repository.
SPATIAL_CONFIG=local_config.yaml make lab
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
