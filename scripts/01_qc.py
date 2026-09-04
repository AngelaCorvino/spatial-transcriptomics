#!/usr/bin/env python3
"""Run QC pipeline using reusable functions from src."""

from __future__ import annotations

import argparse

from spatial_transcriptomics.config import (
    get_config_path,
    get_config_section,
    get_section_path,
    load_config,
)
from spatial_transcriptomics.data import calc_qc_metrics, load_data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run QC on spatial transcriptomics inputs.",
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
    output_dir = get_section_path(config, paths, "qc_output", "results/qc")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not data_path.exists():
        raise FileNotFoundError(f"Data directory not found: {data_path}")

    print("Loading data for QC...")
    adata = load_data(data_path)

    print("Calculating QC metrics...")
    adata = calc_qc_metrics(
        adata,
        mt_genes=config.get("mt_genes"),
        rp_genes=config.get("rp_genes"),
    )

    output_file = output_dir / "qc_data.h5ad"
    print(f"Saving QC output to {output_file}")
    adata.write(output_file)

    print("QC pipeline complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
