"""Analysis helpers for spatial transcriptomics notebooks and scripts."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

# Same inclusive thresholds as the Visium HD notebook; no filtering is applied.
HD_QC_CANDIDATES = (
    ("No additional filter", 0, 0, 100.0),
    ("Lenient 8 µm", 10, 10, 25.0),
    ("Moderate 8 µm", 25, 20, 20.0),
    ("Higher-content 8 µm", 50, 40, 20.0),
)

# Selected on 19 October 2026. Keep v1 fixed for reproducible downstream inputs.
HD_QC_FILTER_VERSION = "thymus_8um_qc_v1"
HD_QC_SELECTED_CANDIDATE = "Moderate 8 µm"


def selected_hd_qc_policy() -> dict[str, Any]:
    """Describe the selected first-pass filter independently of generic defaults."""
    candidate = next(
        row for row in HD_QC_CANDIDATES if row[0] == HD_QC_SELECTED_CANDIDATE
    )
    _, min_counts, min_genes, max_mt = candidate
    if (min_counts, min_genes, max_mt) != (25, 20, 20.0):
        raise ValueError("The v1 QC thresholds have changed; define a new version.")
    return {
        "version": HD_QC_FILTER_VERSION,
        "selected_date": "2026-10-19",
        "candidate": HD_QC_SELECTED_CANDIDATE,
        "bin_size_um": 8,
        "input": "Space Ranger filtered, in-tissue 8 µm bins",
        "min_counts": min_counts,
        "min_genes": min_genes,
        "max_pct_mt": max_mt,
        "thresholds_inclusive": True,
        "anatomical_roi_applied": False,
        "ribosomal_filter_applied": False,
        "upper_counts_or_genes_filter_applied": False,
        "normalization_applied": False,
        "gene_features_filtered": False,
        "scope": "All input tissue; no color, edge, cluster or sample exclusion",
    }


def select_hd_qc_bins(qc: pd.DataFrame) -> pd.DataFrame:
    """Return auditable v1 flags for every input bin, preserving barcode order."""
    policy = selected_hd_qc_policy()
    _, masks = evaluate_hd_qc(qc)
    if not qc.index.is_unique or qc.index.hasnans:
        raise ValueError("Input bin barcodes must be unique and nonmissing.")
    if (qc[["total_counts", "n_genes_by_counts", "pct_counts_mt"]] < 0).any().any():
        raise ValueError("QC metrics cannot be negative.")
    if (qc["pct_counts_mt"] > 100).any():
        raise ValueError("Mitochondrial percentages cannot exceed 100.")
    result = qc.copy()
    result["fails_umi"] = qc["total_counts"] < policy["min_counts"]
    result["fails_genes"] = qc["n_genes_by_counts"] < policy["min_genes"]
    result["fails_mt"] = qc["pct_counts_mt"] > policy["max_pct_mt"]
    result["qc_keep"] = masks[HD_QC_SELECTED_CANDIDATE]
    result["qc_filter_version"] = HD_QC_FILTER_VERSION
    return result


def evaluate_hd_qc(
    qc: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Compare notebook thresholds on all bins, leaving the input unchanged."""
    columns = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing = set(columns) - set(qc.columns)
    if missing:
        raise ValueError(
            f"Missing QC metrics (check reference genes): {sorted(missing)}"
        )
    if qc.empty or not np.isfinite(qc[columns].to_numpy()).all():
        raise ValueError("QC metrics must contain bins and finite values.")
    counts = qc["total_counts"].to_numpy()
    genes = qc["n_genes_by_counts"].to_numpy()
    mt = qc["pct_counts_mt"].to_numpy()
    total = counts.sum()
    rows = []
    masks = {}
    for name, min_counts, min_genes, max_mt in HD_QC_CANDIDATES:
        mask = (counts >= min_counts) & (genes >= min_genes) & (mt <= max_mt)
        masks[name] = mask
        rows.append(
            {
                "candidate": name,
                "min_counts": min_counts,
                "min_genes": min_genes,
                "max_pct_mt": max_mt,
                "bins_retained": int(mask.sum()),
                "bins_retained_pct": 100 * mask.mean(),
                "transcripts_retained_pct": 100 * counts[mask].sum() / total
                if total
                else 0,
                "median_counts_retained": float(np.median(counts[mask]))
                if mask.any()
                else np.nan,
                "median_genes_retained": float(np.median(genes[mask]))
                if mask.any()
                else np.nan,
            }
        )
    return pd.DataFrame(rows), masks


def summarize_hd_qc(qc: pd.DataFrame) -> dict[str, int | float]:
    """Summarize a mouse using the notebook's cross-mouse QC statistics."""
    result: dict[str, int | float] = {"bins": len(qc)}
    for column, label in [
        ("total_counts", "counts"),
        ("n_genes_by_counts", "genes"),
        ("pct_counts_mt", "pct_mt"),
    ]:
        result[f"median_{label}"] = float(qc[column].median())
        for percentile in [95, 99] if label == "pct_mt" else [5, 95]:
            result[f"{label}_p{percentile:02d}"] = float(
                qc[column].quantile(percentile / 100),
            )
    return result


def decompose_hd_qc(qc: pd.DataFrame) -> pd.DataFrame:
    """Count independent and overlapping failures without filtering input bins.

    For each existing candidate, the eight ``exclusive`` groups partition the
    input. The three ``marginal`` groups count all failures of one metric and
    overlap with one another; they must not be added to the exclusive totals.
    Percentages use all input bins or all input UMIs as their denominator.
    """
    evaluate_hd_qc(qc)  # Use the same required-column and finite-value checks.
    counts = qc["total_counts"].to_numpy(dtype=float)
    genes = qc["n_genes_by_counts"].to_numpy()
    mt = qc["pct_counts_mt"].to_numpy()
    total_umis = counts.sum()
    labels = (
        "passes_all",
        "umi_only",
        "genes_only",
        "umi_and_genes_only",
        "mt_only",
        "umi_and_mt_only",
        "genes_and_mt_only",
        "all_three",
    )
    rows = []
    for candidate, min_counts, min_genes, max_mt in HD_QC_CANDIDATES:
        failures = (counts < min_counts, genes < min_genes, mt > max_mt)
        codes = np.zeros(len(qc), dtype=np.uint8)
        for i, mask in enumerate(failures):
            codes[mask] |= 1 << i
        groups = [("exclusive", name, codes == i) for i, name in enumerate(labels)]
        groups += [
            ("marginal", name, mask)
            for name, mask in zip(("fails_umi", "fails_genes", "fails_mt"), failures)
        ]
        for group_type, group, mask in groups:
            umis = float(counts[mask].sum())
            rows.append(
                {
                    "candidate": candidate,
                    "min_counts": min_counts,
                    "min_genes": min_genes,
                    "max_pct_mt": max_mt,
                    "group_type": group_type,
                    "group": group,
                    "bins": int(mask.sum()),
                    "bins_pct": 100 * mask.mean(),
                    "umis": umis,
                    "umis_pct": 100 * umis / total_umis if total_umis else np.nan,
                    "input_bins": len(qc),
                    "input_umis": float(total_umis),
                }
            )
    return pd.DataFrame(rows)


def _load_scanpy() -> Any:
    """Import scanpy lazily so basic package imports stay lightweight."""
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("scanpy is required for analysis helpers.") from exc
    return sc


def cluster_and_umap(
    adata: Any,
    resolution: float = 1.0,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    *,
    n_top_genes: int = 2000,
    exclude_genes: list[str] | None = None,
    random_state: int = 0,
) -> Any:
    """Cluster log-normalized expression; preserve X and count layers.

    Scale only a sparse HVG copy for PCA. Gene exclusions affect the embedding,
    not the stored expression matrix or marker-gene inspection.
    """
    import anndata as ad

    sc = _load_scanpy()
    if adata.n_obs < 4 or min(n_top_genes, n_pcs, n_neighbors) < 2:
        raise ValueError("Clustering requires >=4 bins and dimensions >=2.")
    if resolution <= 0:
        raise ValueError("Clustering resolution must be positive.")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")
    hvg = adata.var["highly_variable"].to_numpy().copy()
    hvg &= ~adata.var_names.isin(exclude_genes or [])
    adata.var["highly_variable"] = hvg
    n_pcs = min(n_pcs, adata.n_obs - 1, int(hvg.sum()) - 1)
    if n_pcs < 2:
        raise ValueError("At least three variable genes are needed for PCA.")
    # Construct only the HVG matrix: do not duplicate raw count layers.
    embedding = ad.AnnData(adata.X[:, hvg].copy())
    sc.pp.scale(embedding, zero_center=False, max_value=10)
    sc.pp.pca(embedding, n_comps=n_pcs, svd_solver="arpack", random_state=random_state)
    adata.obsm["X_pca"] = embedding.obsm["X_pca"]
    adata.varm["PCs"] = np.zeros((adata.n_vars, n_pcs), dtype=np.float32)
    adata.varm["PCs"][hvg] = embedding.varm["PCs"]
    adata.uns["pca"] = embedding.uns["pca"]
    del embedding
    sc.pp.neighbors(
        adata,
        n_neighbors=min(n_neighbors, adata.n_obs - 1),
        n_pcs=n_pcs,
        use_rep="X_pca",
        random_state=random_state,
    )
    cluster_leiden(adata, resolution=resolution, random_state=random_state)
    sc.tl.umap(adata, min_dist=0.1, spread=1.0, random_state=random_state)
    return adata


def cluster_leiden(adata, *, resolution, random_state=0):
    """Recluster the existing neighbor graph without changing its embedding."""
    sc = _load_scanpy()
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="leiden",
        flavor="igraph",
        directed=False,
        n_iterations=2,
        random_state=random_state,
    )


def summarize_annotation_markers(adata: Any, panels: dict[str, list[str]]):
    """Report marker coverage and cluster expression without assigning identities."""
    coverage = []
    expression = []
    clusters = adata.obs["leiden"].cat.categories
    for panel, genes in panels.items():
        for gene in genes:
            present = gene in adata.var_names
            values = np.array([], dtype=float)
            if present:
                matrix = adata[:, gene].X
                values = np.asarray(
                    matrix.toarray() if hasattr(matrix, "toarray") else matrix
                ).ravel()
            coverage.append(
                {
                    "panel": panel,
                    "gene": gene,
                    "in_matrix": present,
                    "bins_expressing": int((values > 0).sum()) if present else np.nan,
                }
            )
            if not present:
                continue
            for cluster in clusters:
                selected = values[adata.obs["leiden"].eq(cluster).to_numpy()]
                expression.append(
                    {
                        "cluster": str(cluster),
                        "panel": panel,
                        "gene": gene,
                        "mean_log_expression": float(selected.mean()),
                        "bins_expressing_pct": float(100 * (selected > 0).mean()),
                    }
                )
    return pd.DataFrame(coverage), pd.DataFrame(expression)


def rank_annotation_markers(adata: Any, n_genes: int = 50) -> pd.DataFrame:
    """Descriptive cluster-vs-rest markers, not treatment-level inference."""
    sc = _load_scanpy()
    sizes = adata.obs["leiden"].value_counts()
    groups = [str(group) for group, size in sizes.items() if 1 < size < adata.n_obs]
    columns = ["group", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    if len(sizes) < 2 or not groups:
        return pd.DataFrame(columns=columns)
    sc.tl.rank_genes_groups(
        adata,
        "leiden",
        groups=groups,
        reference="rest",
        use_raw=False,
        method="t-test_overestim_var",
        n_genes=min(n_genes, adata.n_vars),
        pts=True,
    )
    return sc.get.rank_genes_groups_df(adata, group=None)


def compute_statistics(adata: Any) -> dict[str, int | float]:
    """Compute a small set of summary statistics from an AnnData-like object."""
    stats: dict[str, int | float] = {
        "n_obs": int(getattr(adata, "n_obs", 0)),
        "n_vars": int(getattr(adata, "n_vars", 0)),
    }

    total_counts = getattr(getattr(adata, "obs", {}), "get", lambda *_: None)(
        "total_counts",
    )
    if total_counts is not None:
        stats["mean_counts"] = float(total_counts.mean())
    else:
        stats["mean_counts"] = 0.0

    if hasattr(adata, "obs") and "leiden" in adata.obs:
        stats["n_clusters"] = int(adata.obs["leiden"].nunique())

    return stats


def prepare_results(adata: Any, stats: dict[str, int | float]) -> dict[str, Any]:
    """Prepare the notebook analysis payload for downstream export."""
    return {
        "data": adata,
        "statistics": stats,
        "summary": (
            f"Analysis complete for {stats.get('n_obs', 0)} observations and "
            f"{stats.get('n_vars', 0)} genes."
        ),
    }
