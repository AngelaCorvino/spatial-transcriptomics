#!/usr/bin/env python3
"""Run preprocessing using reusable functions from src."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import re
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

from spatial_transcriptomics.analysis import (
    HD_QC_FILTER_VERSION,
    decompose_hd_qc,
    select_hd_qc_bins,
    selected_hd_qc_policy,
)
from spatial_transcriptomics.config import (
    get_config_path,
    get_config_section,
    get_section_path,
    load_config,
)
from spatial_transcriptomics.data import (
    apply_hd_qc_filter,
    load_data,
    load_reference_genes,
    load_visium_hd_bin,
    preprocess_data,
    stage_visium_hd_qc,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run preprocessing on spatial transcriptomics data.",
    )
    parser.add_argument(
        "--config",
        default="configs/local.yaml",
        help="Path to YAML config file.",
    )
    parser.add_argument(
        "--visium-hd",
        action="store_true",
        help="Apply selected QC v1 to 8 µm bins; save unnormalized retained counts.",
    )
    samples = parser.add_mutually_exclusive_group()
    samples.add_argument("--mice", nargs="+", help="Mouse IDs from FD1 through FD12.")
    samples.add_argument("--all-mice", action="store_true", help="Filter FD1–FD12.")
    parser.add_argument(
        "--source-root", type=Path, help="Parent of FD*/binned_outputs.tar.gz."
    )
    parser.add_argument(
        "--qc-dir", type=Path, help="Parent of saved FD*/bin_qc.parquet."
    )
    parser.add_argument("--scratch-dir", type=Path, help="Default: job TMPDIR.")
    parser.add_argument("--output-dir", type=Path, help="An unused output directory.")
    return parser.parse_args()


def verify_saved_hd_qc(selection: pd.DataFrame, saved_path: Path) -> int:
    """Check barcode identities, metrics, geometry and decisions against the review."""
    saved = pd.read_parquet(saved_path)
    if (
        not saved.index.is_unique
        or len(saved) != len(selection)
        or not selection.index.isin(saved.index).all()
    ):
        raise ValueError(f"Input barcodes differ from the reviewed QC: {saved_path}")
    saved = saved.reindex(selection.index)
    for column in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "x", "y"]:
        if column not in saved or not np.allclose(
            selection[column],
            saved[column],
            rtol=1e-6 if column == "pct_counts_mt" else 0,
            atol=1e-5 if column == "pct_counts_mt" else 0,
        ):
            raise ValueError(f"{column} differs from the reviewed QC: {saved_path}")
    changed = selection["qc_keep"] != select_hd_qc_bins(saved)["qc_keep"]
    # Old float32 ratios could place exactly 20% just above the cutoff.
    rounding_fix = selection["qc_keep"] & selection["pct_counts_mt"].eq(20)
    if (changed & ~rounding_fix).any():
        raise ValueError(
            f"Retained bins differ from the reviewed candidate: {saved_path}"
        )
    return int(changed.sum())


def filter_hd_mouse(mouse, archive, qc_path, output_dir, scratch, mt_genes):
    """Filter one mouse and save its raw counts, all-bin flags and summary."""
    with tempfile.TemporaryDirectory(prefix=f"filter_{mouse}_", dir=scratch) as stage:
        bin_dir = stage_visium_hd_qc(
            archive,
            stage,
            use_filtered_matrix=True,
            include_scalefactors=True,
        )
        adata = load_visium_hd_bin(bin_dir, library_id=mouse, load_images=False)
        adata.obs["sample_id"] = mouse
        retained, selection = apply_hd_qc_filter(adata, mt_genes)
        del adata
    boundary_corrections = verify_saved_hd_qc(selection, qc_path)
    # Refer to the durable archive, not the now-deleted extraction directory.
    for metadata in [
        retained.uns["visium_hd"],
        retained.uns["spatial"][mouse]["metadata"],
    ]:
        metadata["source_path"] = str(archive)
        metadata["archive_member"] = "square_008um/filtered_feature_bc_matrix.h5"
    sample_dir = output_dir / mouse
    sample_dir.mkdir()
    temporary_h5ad = sample_dir / "filtered_counts.partial.h5ad"
    retained.write_h5ad(temporary_h5ad, compression="gzip")
    temporary_h5ad.replace(sample_dir / "filtered_counts.h5ad")
    selection.to_parquet(sample_dir / "bin_selection.parquet")
    keep = selection["qc_keep"]
    input_umis = float(selection["total_counts"].sum())
    retained_umis = float(selection.loc[keep, "total_counts"].sum())
    failure = decompose_hd_qc(selection)
    failure = failure[failure["candidate"] == selected_hd_qc_policy()["candidate"]]
    failure.to_csv(sample_dir / "failure_summary.csv", index=False)
    return {
        "mouse_id": mouse,
        "input_bins": len(selection),
        "bins_retained": int(keep.sum()),
        "bins_excluded": int((~keep).sum()),
        "bins_excluded_pct": float(100 * (~keep).mean()),
        "input_umis": input_umis,
        "umis_retained": retained_umis,
        "umis_excluded_pct": 100 * (input_umis - retained_umis) / input_umis,
        "features_retained": retained.n_vars,
        "mt_genes_matched": len(retained.uns["qc_filter"]["mt_genes_matched"]),
        "mt_boundary_corrections": boundary_corrections,
        "source_archive": str(archive),
        "reviewed_qc": str(qc_path),
        "reviewed_qc_sha256": hashlib.sha256(qc_path.read_bytes()).hexdigest(),
    }


def run_visium_hd_filter(config: dict[str, object], args: argparse.Namespace) -> int:
    """Apply QC v1 sequentially on SLURM, preserving existing data and outputs."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Apply experimental-data filtering in a SLURM compute job.")
    hd = get_config_section(config, "visium_hd")
    if int(str(hd.get("bin_size_um", 8))) != 8 or not hd.get(
        "use_filtered_matrix", True
    ):
        raise ValueError("QC v1 requires Space Ranger filtered 8 µm inputs.")
    mice = [f"FD{i}" for i in range(1, 13)] if args.all_mice else args.mice
    if (
        not mice
        or len(mice) != len(set(mice))
        or any(re.fullmatch(r"FD(?:[1-9]|1[0-2])", mouse) is None for mouse in mice)
    ):
        raise ValueError("Specify --all-mice or unique --mice FD1 … FD12.")
    paths = get_config_section(config, "paths")
    repo_root = get_config_path(config, "repo_root", ".")
    source_root = (
        args.source_root
        or get_section_path(config, hd, "source_dir", "data/FD1").parent
    )
    qc_dir = (
        args.qc_dir
        or get_section_path(config, paths, "qc_output", "results/qc")
        / "visium_hd_008um"
    )
    output_dir = (
        args.output_dir
        or get_section_path(config, paths, "preprocess_output", "results/preprocess")
        / "visium_hd_008um"
        / HD_QC_FILTER_VERSION
    )
    source_root, qc_dir, output_dir = [
        p.expanduser().resolve() for p in (source_root, qc_dir, output_dir)
    ]
    for protected in [repo_root, source_root, qc_dir]:
        if output_dir == protected or protected in output_dir.parents:
            raise ValueError(
                f"Keep filter outputs outside code, raw data and QC: {protected}"
            )
    if output_dir.exists():
        raise FileExistsError(
            f"Preserve existing results; choose an unused --output-dir: {output_dir}"
        )
    scratch_value = args.scratch_dir or os.environ.get("TMPDIR")
    if not scratch_value:
        raise ValueError("Set TMPDIR in the compute job or supply --scratch-dir.")
    scratch = Path(scratch_value).expanduser().resolve()
    mt_path = repo_root / "reference_genomes" / "mouse_mitochondrial_genes.txt"
    mt_genes = load_reference_genes(mt_path)
    if not mt_genes:
        raise ValueError(f"Mitochondrial reference list missing or empty: {mt_path}")
    archives = {mouse: source_root / mouse / "binned_outputs.tar.gz" for mouse in mice}
    for mouse, archive in archives.items():
        for required in [archive, qc_dir / mouse / "bin_qc.parquet"]:
            if not required.is_file():
                raise FileNotFoundError(required)
    policy = selected_hd_qc_policy()
    run = {
        "status": "running",
        "filter": policy,
        "mice": mice,
        "completed_mice": [],
        "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source_root": str(source_root),
        "qc_dir": str(qc_dir),
        "output_dir": str(output_dir),
        "resolved_config": config,
        "source_sha256": {
            relative: hashlib.sha256((repo_root / relative).read_bytes()).hexdigest()
            for relative in [
                "scripts/02_preprocess.py",
                "src/spatial_transcriptomics/analysis.py",
                "src/spatial_transcriptomics/data.py",
                "reference_genomes/mouse_mitochondrial_genes.txt",
            ]
        },
    }
    scratch.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True)
    run_path = output_dir / "run.json"
    run_path.write_text(json.dumps(run, indent=2) + "\n")
    summaries = []
    started = time.monotonic()
    try:
        for mouse, archive in archives.items():
            print(f"{mouse}: applying {HD_QC_FILTER_VERSION}", flush=True)
            summary = filter_hd_mouse(
                mouse,
                archive,
                qc_dir / mouse / "bin_qc.parquet",
                output_dir,
                scratch,
                mt_genes,
            )
            summaries.append(summary)
            pd.DataFrame(summaries).to_csv(
                output_dir / "filter_summary.csv", index=False
            )
            run["completed_mice"].append(mouse)
            run_path.write_text(json.dumps(run, indent=2) + "\n")
            print(
                f"{mouse}: retained {summary['bins_retained']:,} "
                f"/ {summary['input_bins']:,} bins",
                flush=True,
            )
            gc.collect()
        run["status"] = "complete"
    except Exception as exc:
        run.update(status="failed", error=str(exc))
        raise
    finally:
        run["elapsed_seconds"] = time.monotonic() - started
        run_path.write_text(json.dumps(run, indent=2) + "\n")
    print(f"QC filter complete. Retained raw counts: {output_dir}", flush=True)
    return 0


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    if args.visium_hd:
        return run_visium_hd_filter(config, args)

    data_path = get_config_path(config, "data_dir", ".")
    paths = get_config_section(config, "paths")
    output_dir = args.output_dir or get_section_path(
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
