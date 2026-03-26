# Reference Genome Lists

This folder contains gene lists and reference annotations for filtering and analysis.

## Contents

### `mouse_mitochondrial_genes.txt`
List of mouse mitochondrial gene symbols (mt-*).
Used to calculate mitochondrial QC metrics.

### `mouse_ribosomal_genes.txt`
List of mouse ribosomal protein genes (Rpl*, Rps*).
Used to identify ribosomal contamination in spatial data.

### `mouse_protein_coding.txt`
List of mouse protein-coding genes (Ensembl).
Use for filtering to protein-coding genes only.

## Usage Example

```python
import pandas as pd

# Load mitochondrial genes
mt_genes = pd.read_csv('reference_genomes/mouse_mitochondrial_genes.txt', header=None)[0].tolist()

# Filter to genes in data
mt_genes_in_data = [g for g in mt_genes if g in adata.var_names]

# Calculate fraction of counts from mt genes
adata.obs['pct_counts_mt'] = (
    adata[:, mt_genes_in_data].X.sum(axis=1) / adata.X.sum(axis=1) * 100
)
```

## Sources

- **Mitochondrial genes:** NCBI Gene database, prefix `mt-`
- **Ribosomal proteins:** GO term GO:0022625 (proteasome complex)
- **Protein-coding:** Ensembl database

## Updating Gene Lists

When updating reference genomes (e.g., for a new mouse build or species), ensure:
1. One gene symbol per line
2. Use official gene symbols (e.g., `AC3`, `ACTIN`)
3. Case-sensitive (typically lowercase for mouse)
4. Remove duplicates
