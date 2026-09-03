"""Data loading and preprocessing helpers for spatial transcriptomics."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _load_scanpy() -> Any:
    """Import scanpy lazily so the package can be imported without it."""
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError(
            "scanpy is required for data loading and preprocessing helpers.",
        ) from exc
    return sc


def _load_squidpy() -> Any:
    """Import squidpy lazily so notebook-specific workflows stay optional."""
    try:
        import squidpy as sq
    except ImportError as exc:
        raise ImportError(
            "squidpy is required to read Visium directories.",
        ) from exc
    return sq


def load_reference_genes(gene_file: str | Path) -> list[str]:
    """Load a newline-delimited gene list from disk."""
    path = Path(gene_file)
    if not path.exists():
        return []

    with path.open() as handle:
        return [line.strip() for line in handle if line.strip()]


def load_data(data_path: str | Path) -> Any:
    """Load spatial transcriptomics data from a file or Visium directory.

    Supported inputs:
    - `.h5ad` files via `scanpy.read_h5ad`
    - Visium output directories via `squidpy.read.visium`
    """
    path = Path(data_path)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.is_file() and path.suffix == ".h5ad":
        sc = _load_scanpy()
        adata = sc.read_h5ad(path)
    elif path.is_dir():
        sq = _load_squidpy()
        adata = sq.read.visium(path)
    else:
        message = (
            f"Unsupported data path: {path}. Expected a .h5ad file or Visium folder."
        )
        raise ValueError(
            message,
        )

    adata.var_names_make_unique()
    return adata


def load_visium_hd_bin(
    data_path: str | Path,
    bin_size_um: int = 8,
    library_id: str | None = None,
    load_images: bool = True,
    use_filtered_matrix: bool = True,
) -> Any:
    """Load one extracted Visium HD bin size into an AnnData object.

    ``data_path`` may be the staged Space Ranger root, its ``binned_outputs``
    directory, or the selected ``square_XXXum`` directory itself. The function
    attaches Space Ranger's bin metadata and full-resolution pixel coordinates
    without applying any downstream QC filters or normalization.
    """
    if bin_size_um <= 0:
        raise ValueError("bin_size_um must be a positive integer.")

    path = Path(data_path).expanduser()
    bin_name = f"square_{bin_size_um:03d}um"
    candidates = [
        path,
        path / bin_name,
        path / "binned_outputs" / bin_name,
    ]
    bin_path = next(
        (
            candidate
            for candidate in candidates
            if candidate.name == bin_name and candidate.is_dir()
        ),
        None,
    )
    if bin_path is None:
        expected = path / "binned_outputs" / bin_name
        raise FileNotFoundError(
            f"Visium HD {bin_size_um} µm directory not found. Expected {expected}",
        )

    matrix_name = (
        "filtered_feature_bc_matrix.h5"
        if use_filtered_matrix
        else "raw_feature_bc_matrix.h5"
    )
    matrix_path = bin_path / matrix_name
    spatial_path = bin_path / "spatial"
    positions_path = spatial_path / "tissue_positions.parquet"

    if not matrix_path.is_file():
        raise FileNotFoundError(f"Visium HD count matrix not found: {matrix_path}")
    if not positions_path.is_file():
        raise FileNotFoundError(
            f"Visium HD tissue positions not found: {positions_path}",
        )

    sc = _load_scanpy()
    adata = sc.read_10x_h5(matrix_path, gex_only=True)
    adata.var_names_make_unique()
    adata.obs_names = adata.obs_names.astype(str)

    positions = pd.read_parquet(positions_path)
    if "barcode" not in positions.columns:
        raise ValueError(f"{positions_path} does not contain a 'barcode' column.")
    if positions["barcode"].duplicated().any():
        raise ValueError(f"{positions_path} contains duplicate barcodes.")

    positions = positions.copy()
    positions["barcode"] = positions["barcode"].astype(str)
    positions = positions.set_index("barcode").reindex(adata.obs_names)

    x_column = "pxl_col_in_fullres"
    y_column = "pxl_row_in_fullres"
    missing_columns = [
        column for column in (x_column, y_column) if column not in positions.columns
    ]
    if missing_columns:
        raise ValueError(
            f"{positions_path} is missing coordinate columns {missing_columns}.",
        )
    if positions[[x_column, y_column]].isna().to_numpy().any():
        raise ValueError("Some matrix barcodes do not have spatial coordinates.")

    for column in positions.columns:
        adata.obs[column] = positions[column].to_numpy()
    adata.obs["bin_size_um"] = bin_size_um
    adata.obsm["spatial"] = positions[[x_column, y_column]].to_numpy(dtype=float)

    scalefactors_path = spatial_path / "scalefactors_json.json"
    scalefactors: dict[str, object] = {}
    if scalefactors_path.is_file():
        with scalefactors_path.open() as handle:
            scalefactors = json.load(handle)

    images: dict[str, Any] = {}
    if load_images:
        try:
            import matplotlib.image as mpimg
        except ImportError as exc:
            raise ImportError("matplotlib is required to load Visium images.") from exc

        for image_key, image_name in {
            "hires": "tissue_hires_image.png",
            "lowres": "tissue_lowres_image.png",
        }.items():
            image_path = spatial_path / image_name
            if image_path.is_file():
                images[image_key] = mpimg.imread(image_path)

    resolved_library_id = library_id or f"{path.name}_{bin_name}"
    adata.uns.setdefault("spatial", {})[resolved_library_id] = {
        "images": images,
        "scalefactors": scalefactors,
        "metadata": {
            "bin_size_um": bin_size_um,
            "matrix": matrix_name,
            "source_path": str(bin_path.resolve()),
        },
    }
    adata.uns["visium_hd"] = {
        "bin_size_um": bin_size_um,
        "library_id": resolved_library_id,
        "source_path": str(bin_path.resolve()),
    }
    return adata


def calc_qc_metrics(
    adata: Any,
    mt_genes: list[str] | None = None,
    rp_genes: list[str] | None = None,
) -> Any:
    """Calculate standard QC metrics plus MT/RP percentages when provided."""
    sc = _load_scanpy()
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    total_counts = np.asarray(adata.X.sum(axis=1)).ravel()

    if mt_genes:
        mt_genes_in_data = [gene for gene in mt_genes if gene in adata.var_names]
        if mt_genes_in_data:
            mt_counts = np.asarray(adata[:, mt_genes_in_data].X.sum(axis=1)).ravel()
            adata.obs["pct_counts_mt"] = np.divide(
                mt_counts,
                total_counts,
                out=np.zeros_like(mt_counts, dtype=float),
                where=total_counts > 0,
            ) * 100

    if rp_genes:
        rp_genes_in_data = [gene for gene in rp_genes if gene in adata.var_names]
        if rp_genes_in_data:
            rp_counts = np.asarray(adata[:, rp_genes_in_data].X.sum(axis=1)).ravel()
            adata.obs["pct_counts_rp"] = np.divide(
                rp_counts,
                total_counts,
                out=np.zeros_like(rp_counts, dtype=float),
                where=total_counts > 0,
            ) * 100

    return adata


def filter_by_quality(
    adata: Any,
    min_counts: int = 500,
    max_counts: int | None = None,
    min_genes: int = 250,
    max_genes: int | None = None,
    max_pct_mt: float = 20.0,
    max_pct_rp: float = 30.0,
) -> tuple[Any, pd.DataFrame]:
    """Filter spots by QC metrics and return the filtered object plus summary."""
    sc = _load_scanpy()

    if "total_counts" not in adata.obs or "n_genes_by_counts" not in adata.obs:
        adata = calc_qc_metrics(adata)

    before = {
        "n_obs": adata.n_obs,
        "n_counts": float(adata.obs["total_counts"].sum()),
        "mean_counts": float(adata.obs["total_counts"].mean()),
    }

    filtered = adata.copy()
    sc.pp.filter_cells(filtered, min_counts=min_counts, inplace=True)
    if max_counts is not None:
        sc.pp.filter_cells(filtered, max_counts=max_counts, inplace=True)

    sc.pp.filter_cells(filtered, min_genes=min_genes, inplace=True)
    if max_genes is not None:
        sc.pp.filter_cells(filtered, max_genes=max_genes, inplace=True)

    if "pct_counts_mt" in filtered.obs:
        filtered = filtered[filtered.obs["pct_counts_mt"] < max_pct_mt].copy()

    if "pct_counts_rp" in filtered.obs:
        filtered = filtered[filtered.obs["pct_counts_rp"] < max_pct_rp].copy()

    after = {
        "n_obs": filtered.n_obs,
        "n_counts": float(filtered.obs["total_counts"].sum()),
        "mean_counts": float(filtered.obs["total_counts"].mean()),
    }

    summary = pd.DataFrame([before, after], index=["Before", "After"])
    return filtered, summary


def normalize_and_log(adata: Any, target_sum: float = 1e4) -> Any:
    """Normalize counts, preserve raw counts, and apply log1p base 2."""
    sc = _load_scanpy()

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)
    sc.pp.log1p(adata, base=2, inplace=True)
    return adata


def preprocess_data(
    adata: Any,
    mt_genes: list[str] | None = None,
    rp_genes: list[str] | None = None,
    min_counts: int = 500,
    max_counts: int | None = None,
    min_genes: int = 250,
    max_genes: int | None = None,
    max_pct_mt: float = 20.0,
    max_pct_rp: float = 30.0,
    target_sum: float = 1e4,
) -> Any:
    """Run the standard preprocessing pipeline used in the notebook."""
    adata = calc_qc_metrics(adata, mt_genes=mt_genes, rp_genes=rp_genes)
    adata, _ = filter_by_quality(
        adata,
        min_counts=min_counts,
        max_counts=max_counts,
        min_genes=min_genes,
        max_genes=max_genes,
        max_pct_mt=max_pct_mt,
        max_pct_rp=max_pct_rp,
    )
    return normalize_and_log(adata, target_sum=target_sum)


def read_de_csv(path: str | Path) -> pd.DataFrame:
    """Read DESeq2 results CSV."""
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)

    required = {"stat", "padj"}
    missing = required.difference(df.columns)
    if missing:
        message = f"{Path(path).name}: missing columns {sorted(missing)}"
        raise ValueError(message)

    return (
        df.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["stat", "padj"])
        .loc[lambda data: ~data.index.duplicated(keep="first")]
    )
