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
) -> Any:
    """Run HVG selection, PCA, neighborhood graph, Leiden, and UMAP."""
    sc = _load_scanpy()

    sc.pp.highly_variable_genes(adata, inplace=True)
    sc.pp.scale(adata, inplace=True)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X_pca")
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden", inplace=True)
    sc.tl.umap(adata, min_dist=0.1, spread=1.0)
    return adata


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
