"""Spatial transcriptomics package."""

from .analysis import cluster_and_umap, compute_statistics, prepare_results
from .config import get_data_path, get_output_path, load_config
from .data import (
    calc_qc_metrics,
    filter_by_quality,
    load_data,
    load_reference_genes,
    normalize_and_log,
    preprocess_data,
    read_de_csv,
)
from .plotting import pretty_title, stacked_bar, summary_dotplot

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "calc_qc_metrics",
    "cluster_and_umap",
    "compute_statistics",
    "filter_by_quality",
    "get_data_path",
    "get_output_path",
    "load_config",
    "load_data",
    "load_reference_genes",
    "normalize_and_log",
    "pretty_title",
    "prepare_results",
    "preprocess_data",
    "read_de_csv",
    "stacked_bar",
    "summary_dotplot",
]
