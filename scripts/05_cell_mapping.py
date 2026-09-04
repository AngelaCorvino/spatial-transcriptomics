#!/usr/bin/env python3
"""Run cell mapping workflow and save outputs."""

from __future__ import annotations

import argparse

from spatial_transcriptomics.config import (
    get_config_section,
    get_section_path,
    load_config,
)
from spatial_transcriptomics.data import load_data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run cell mapping workflow.")
    parser.add_argument(
        "--config",
        default="configs/local.yaml",
        help="Path to YAML config file.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)

    paths = get_config_section(config, "paths")
    input_path = (
        get_section_path(config, paths, "spatial_output", "results/spatial")
        / "spatial_data.h5ad"
    )
    output_dir = get_section_path(
        config,
        paths,
        "cell_mapping_output",
        "results/cell_mapping",
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Spatial analysis input not found: {input_path}")

    print("Loading spatial analysis results for cell mapping...")
    adata = load_data(input_path)

    print("Saving cell mapping placeholder output...")
    output_file = output_dir / "cell_mapping.h5ad"
    adata.write(output_file)

    print("Cell mapping complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
