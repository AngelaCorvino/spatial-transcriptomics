#!/usr/bin/env python3
"""Run data integration using reusable functions from src."""

from __future__ import annotations

import argparse

from spatial_transcriptomics.analysis import (
    cluster_and_umap,
    compute_statistics,
    prepare_results,
)
from spatial_transcriptomics.config import (
    get_config_section,
    get_section_path,
    load_config,
)
from spatial_transcriptomics.data import load_data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run integration and clustering.")
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
        get_section_path(config, paths, "preprocess_output", "results/preprocess")
        / "preprocessed_data.h5ad"
    )
    output_dir = get_section_path(
        config,
        paths,
        "integration_output",
        "results/integration",
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Preprocessed input not found: {input_path}")

    print("Loading preprocessed data for integration...")
    adata = load_data(input_path)

    print("Running clustering and UMAP...")
    analysis = config.get("analysis", {})
    adata = cluster_and_umap(
        adata,
        resolution=analysis.get("resolution", 1.0),
        n_pcs=analysis.get("n_pcs", 30),
        n_neighbors=analysis.get("n_neighbors", 15),
    )

    stats = compute_statistics(adata)
    results = prepare_results(adata, stats)

    output_file = output_dir / "integration_result.h5ad"
    print(f"Saving integration results to {output_file}")
    results["data"].write(output_file)

    summary_file = output_dir / "integration_summary.txt"
    print(f"Writing integration summary to {summary_file}")
    summary_file.write_text(results["summary"] + "\n")

    print("Integration complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
