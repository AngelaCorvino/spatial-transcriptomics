"""Analysis helpers for spatial transcriptomics notebooks and scripts."""

from typing import Any


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


def compute_statistics(adata: Any) -> dict[str, float]:
    """Compute a small set of summary statistics from an AnnData-like object."""
    stats = {
        "n_obs": int(getattr(adata, "n_obs", 0)),
        "n_vars": int(getattr(adata, "n_vars", 0)),
    }

    total_counts = getattr(getattr(adata, "obs", {}), "get", lambda *_: None)(
        "total_counts"
    )
    if total_counts is not None:
        stats["mean_counts"] = float(total_counts.mean())
    else:
        stats["mean_counts"] = 0.0

    if hasattr(adata, "obs") and "leiden" in adata.obs:
        stats["n_clusters"] = int(adata.obs["leiden"].nunique())

    return stats


def prepare_results(adata: Any, stats: dict[str, float]) -> dict[str, Any]:
    """Prepare the notebook analysis payload for downstream export."""
    return {
        "data": adata,
        "statistics": stats,
        "summary": (
            f"Analysis complete for {stats.get('n_obs', 0)} observations and "
            f"{stats.get('n_vars', 0)} genes."
        ),
    }
