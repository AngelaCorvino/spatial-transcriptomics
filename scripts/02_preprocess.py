#!/usr/bin/env python3
"""Run preprocessing using reusable functions from src."""

from __future__ import annotations

import argparse

from spatial_transcriptomics.config import (
    get_config_path,
    get_config_section,
    get_section_path,
    load_config,
)
from spatial_transcriptomics.data import load_data, preprocess_data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run preprocessing on spatial transcriptomics data.",
    )
    parser.add_argument(
        "--config",
        default="configs/local.yaml",
        help="Path to YAML config file.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)

    data_path = get_config_path(config, "data_dir", ".")
    paths = get_config_section(config, "paths")
    output_dir = get_section_path(
        config,
        paths,
        "preprocess_output",
        "results/preprocess",
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if not data_path.exists():
        raise FileNotFoundError(f"Data directory not found: {data_path}")

    print("Loading data for preprocessing...")
    adata = load_data(data_path)

    print("Running preprocessing pipeline...")
    analysis = config.get("analysis", {})
    adata = preprocess_data(
        adata,
        mt_genes=config.get("mt_genes"),
        rp_genes=config.get("rp_genes"),
        min_counts=analysis.get("min_counts", 500),
        max_counts=analysis.get("max_counts"),
        min_genes=analysis.get("min_genes", 250),
        max_genes=analysis.get("max_genes"),
        max_pct_mt=analysis.get("max_pct_mt", 20.0),
        max_pct_rp=analysis.get("max_pct_rp", 30.0),
        target_sum=analysis.get("target_sum", 10000.0),
    )

    output_file = output_dir / "preprocessed_data.h5ad"
    print(f"Saving preprocessed output to {output_file}")
    adata.write(output_file)

    print("Preprocessing complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
