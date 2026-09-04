#!/usr/bin/env python3
"""Run spatial analysis and save results."""

from __future__ import annotations

import argparse

from spatial_transcriptomics.config import (
    get_config_section,
    get_section_path,
    load_config,
)
from spatial_transcriptomics.data import load_data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run spatial analysis workflow.")
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
        get_section_path(config, paths, "integration_output", "results/integration")
        / "integration_result.h5ad"
    )
    output_dir = get_section_path(
        config,
        paths,
        "spatial_output",
        "results/spatial",
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Integration input not found: {input_path}")

    print("Loading integration result for spatial analysis...")
    adata = load_data(input_path)

    print("Saving spatial analysis placeholder output...")
    output_file = output_dir / "spatial_data.h5ad"
    adata.write(output_file)

    print("Spatial analysis complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
