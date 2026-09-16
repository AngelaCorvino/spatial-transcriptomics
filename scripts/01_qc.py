#!/usr/bin/env python3
"""Run QC pipeline using reusable functions from src."""

from __future__ import annotations

import argparse
import gc
import json
import os
import re
import tempfile
import time
from pathlib import Path

import pandas as pd

from spatial_transcriptomics.analysis import evaluate_hd_qc, summarize_hd_qc
from spatial_transcriptomics.config import (
    get_config_path,
    get_config_section,
    get_section_path,
    load_config,
)
from spatial_transcriptomics.data import (
    calc_qc_metrics,
    load_data,
    load_reference_genes,
    load_visium_hd_bin,
    stage_visium_hd_qc,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run QC on spatial transcriptomics inputs.",
    )
    parser.add_argument(
        "--config",
        default="configs/local.yaml",
        help="Path to YAML config file.",
    )
    parser.add_argument(
        "--visium-hd",
        action="store_true",
        help="Run the notebook's 8 µm cross-mouse QC analysis.",
    )
    samples = parser.add_mutually_exclusive_group()
    samples.add_argument("--mice", nargs="+", default=["FD1", "FD2"])
    samples.add_argument("--all-mice", action="store_true", help="Evaluate FD1–FD12.")
    parser.add_argument(
        "--source-root",
        type=Path,
        help="Parent of the FD*/binned_outputs.tar.gz archives.",
    )
    parser.add_argument("--scratch-dir", type=Path, help="Default: job TMPDIR.")
    parser.add_argument(
        "--output-dir", type=Path, help="Override the QC output folder."
    )
    parser.add_argument("--force", action="store_true", help="Recompute cached HD QC.")
    parser.add_argument("--no-plots", action="store_true", help="Write HD tables only.")
    return parser.parse_args()


def run_visium_hd_qc(config: dict[str, object], args: argparse.Namespace) -> int:
    """Run resumable, sequential HD QC without retaining expression matrices."""
    import matplotlib

    matplotlib.use("Agg")
    from spatial_transcriptomics.plotting import plot_hd_qc, plot_hd_qc_comparison

    hd = get_config_section(config, "visium_hd")
    if int(str(hd.get("bin_size_um", 8))) != 8:
        raise ValueError(
            "These notebook candidate thresholds are defined for 8 µm bins."
        )
    mice = [f"FD{i}" for i in range(1, 13)] if args.all_mice else args.mice
    if len(mice) != len(set(mice)) or any(
        re.fullmatch(r"FD(?:[1-9]|1[0-2])", mouse) is None for mouse in mice
    ):
        raise ValueError("Supply unique mouse IDs from FD1 through FD12.")
    source_root = (
        args.source_root
        or get_section_path(
            config,
            hd,
            "source_dir",
            "data/FD1",
        ).parent
    )
    source_root = source_root.expanduser().resolve()
    paths = get_config_section(config, "paths")
    output_dir = args.output_dir or (
        get_section_path(config, paths, "qc_output", "results/qc") / "visium_hd_008um"
    )
    output_dir = output_dir.expanduser().resolve()
    scratch = args.scratch_dir or os.environ.get("TMPDIR")
    if not scratch:
        raise ValueError("Set TMPDIR in a compute job or provide --scratch-dir.")
    scratch = Path(scratch).expanduser().resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    repo_root = get_config_path(config, "repo_root", ".")
    mt_genes = load_reference_genes(
        repo_root / "reference_genomes" / "mouse_mitochondrial_genes.txt",
    )
    rp_genes = load_reference_genes(
        repo_root / "reference_genomes" / "mouse_ribosomal_genes.txt",
    )
    if not mt_genes or not rp_genes:
        raise ValueError(
            "Mouse mitochondrial and ribosomal reference lists are required."
        )
    use_filtered = bool(hd.get("use_filtered_matrix", True))
    archives = {mouse: source_root / mouse / "binned_outputs.tar.gz" for mouse in mice}
    for archive in archives.values():
        if not archive.is_file():
            raise FileNotFoundError(archive)

    comparison_dir = output_dir / "comparisons" / "_".join(mice)
    comparison_dir.mkdir(parents=True, exist_ok=True)
    run_record = {
        "status": "running",
        "mice": mice,
        "config": str(args.config),
        "source_root": str(source_root),
        "bin_size_um": 8,
        "use_filtered_matrix": use_filtered,
        "plots_requested": not args.no_plots,
    }
    run_path = comparison_dir / "run.json"
    run_path.write_text(json.dumps(run_record, indent=2) + "\n")
    started = time.monotonic()
    summaries = []
    candidates = []
    qc_paths = {}
    for mouse, archive in archives.items():
        sample_start = time.monotonic()
        sample_dir = output_dir / mouse
        sample_dir.mkdir(parents=True, exist_ok=True)
        qc_path = sample_dir / "bin_qc.parquet"
        metadata_path = sample_dir / "qc_metadata.json"
        stat = archive.stat()
        fingerprint = {
            "qc_cache_version": 1,
            "archive": str(archive),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "bin_size_um": 8,
            "use_filtered_matrix": use_filtered,
            "mt_genes": mt_genes,
            "rp_genes": rp_genes,
        }
        metadata = {}
        if not args.force and metadata_path.exists():
            try:
                saved_metadata = json.loads(metadata_path.read_text())
            except json.JSONDecodeError:
                print(f"{mouse}: incomplete cache metadata; recomputing QC", flush=True)
            else:
                if isinstance(saved_metadata, dict):
                    metadata = saved_metadata
        if not args.force and qc_path.exists() and metadata.get("input") == fingerprint:
            print(f"{mouse}: reusing cached QC metrics", flush=True)
            qc = pd.read_parquet(qc_path)
        else:
            print(
                f"{mouse}: staging 8 µm matrix and positions from {archive}", flush=True
            )
            with tempfile.TemporaryDirectory(
                prefix=f"qc_{mouse}_", dir=scratch
            ) as stage:
                bin_dir = stage_visium_hd_qc(
                    archive, stage, use_filtered_matrix=use_filtered
                )
                print(f"{mouse}: loading sparse counts and calculating QC", flush=True)
                adata = load_visium_hd_bin(
                    bin_dir,
                    library_id=mouse,
                    load_images=False,
                    use_filtered_matrix=use_filtered,
                )
                mt_matches = sum(gene in adata.var_names for gene in mt_genes)
                rp_matches = sum(gene in adata.var_names for gene in rp_genes)
                if not mt_matches:
                    raise ValueError(
                        f"{mouse}: no mitochondrial reference genes matched."
                    )
                if not rp_matches:
                    print(
                        f"WARNING: {mouse}: no ribosomal reference genes matched.",
                        flush=True,
                    )
                adata = calc_qc_metrics(
                    adata,
                    mt_genes=mt_genes,
                    rp_genes=rp_genes,
                    percent_top=None,
                )
                columns = [
                    column
                    for column in [
                        "total_counts",
                        "n_genes_by_counts",
                        "pct_counts_mt",
                        "pct_counts_rp",
                    ]
                    if column in adata.obs
                ]
                qc = adata.obs[columns].copy()
                qc.index.name = "barcode"
                qc[["x", "y"]] = adata.obsm["spatial"]
                metadata = {
                    "input": fingerprint,
                    "mt_genes_matched": mt_matches,
                    "rp_genes_matched": rp_matches,
                    "features": adata.n_vars,
                }
                del adata
                gc.collect()
            # Validate before publishing a cache; the marker is written last.
            evaluate_hd_qc(qc)
            metadata_path.unlink(missing_ok=True)
            temporary_qc = sample_dir / "bin_qc.tmp.parquet"
            qc.to_parquet(temporary_qc)
            temporary_qc.replace(qc_path)
            temporary_metadata = sample_dir / "qc_metadata.tmp.json"
            temporary_metadata.write_text(json.dumps(metadata, indent=2) + "\n")
            temporary_metadata.replace(metadata_path)

        candidate_summary, _ = evaluate_hd_qc(qc)
        candidate_summary.insert(0, "mouse_id", mouse)
        candidate_summary.to_csv(sample_dir / "candidate_summary.csv", index=False)
        metric_columns = [column for column in qc if column not in {"x", "y"}]
        qc[metric_columns].describe(
            percentiles=[0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99],
        ).T.to_csv(sample_dir / "qc_summary.csv")
        summaries.append(
            {
                "mouse_id": mouse,
                **summarize_hd_qc(qc),
                "mt_genes_matched": metadata["mt_genes_matched"],
                "rp_genes_matched": metadata["rp_genes_matched"],
            }
        )
        candidates.append(candidate_summary)
        qc_paths[mouse] = qc_path
        if not args.no_plots:
            print(f"{mouse}: saving plots", flush=True)
            plot_hd_qc(qc, mouse, sample_dir)
        del qc
        gc.collect()
        print(
            f"{mouse}: finished in {time.monotonic() - sample_start:.1f}s", flush=True
        )

    combined = pd.concat(candidates, ignore_index=True)
    pd.DataFrame(summaries).to_csv(
        comparison_dir / "cross_mouse_qc_summary.csv", index=False
    )
    combined.to_csv(comparison_dir / "cross_mouse_candidate_summary.csv", index=False)
    if not args.no_plots:
        plot_hd_qc_comparison(combined, qc_paths, comparison_dir)
    run_record.update(status="complete", elapsed_seconds=time.monotonic() - started)
    run_path.write_text(json.dumps(run_record, indent=2) + "\n")
    print(f"QC complete. Comparison outputs: {comparison_dir}", flush=True)
    return 0


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    if args.visium_hd:
        return run_visium_hd_qc(config, args)

    data_path = get_config_path(config, "data_dir", ".")
    paths = get_config_section(config, "paths")
    output_dir = args.output_dir or get_section_path(
        config, paths, "qc_output", "results/qc"
    )
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
