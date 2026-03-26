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
from .example import add
from .plotting import pretty_title, stacked_bar, summary_dotplot

__version__ = "0.1.0"

__all__ = [
    "add",
    "load_config",
    "get_data_path",
    "get_output_path",
    "load_data",
    "load_reference_genes",
    "calc_qc_metrics",
    "preprocess_data",
    "filter_by_quality",
    "normalize_and_log",
    "read_de_csv",
    "cluster_and_umap",
    "compute_statistics",
    "prepare_results",
    "stacked_bar",
    "summary_dotplot",
    "pretty_title",
]
