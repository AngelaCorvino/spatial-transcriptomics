#!/usr/bin/env python3
"""Run clustering, including per-mouse Visium HD tissue-annotation review."""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from spatial_transcriptomics.analysis import (
    HD_QC_FILTER_VERSION,
    cluster_and_umap,
    compute_statistics,
    prepare_results,
    rank_annotation_markers,
    select_hd_qc_bins,
    summarize_annotation_markers,
)
from spatial_transcriptomics.config import (
    get_config_path,
    get_config_section,
    get_section_path,
    load_config,
)
from spatial_transcriptomics.data import (
    load_data,
    load_visium_hd_histology,
    normalize_and_log,
)
from spatial_transcriptomics.plotting import plot_hd_annotation, plot_hd_qc_histology


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run integration and clustering.")
    parser.add_argument(
        "--config",
        default="configs/local.yaml",
        help="Path to YAML config file.",
    )
    parser.add_argument(
        "--visium-hd",
        action="store_true",
        help="Cluster filtered 8 µm bins separately for each mouse.",
    )
    samples = parser.add_mutually_exclusive_group()
    samples.add_argument("--mice", nargs="+")
    samples.add_argument("--all-mice", action="store_true")
    parser.add_argument(
        "--input-dir", type=Path, help="Parent of FD*/filtered_counts.h5ad."
    )
    parser.add_argument("--output-dir", type=Path, help="An unused output directory.")
    parser.add_argument("--markers", type=Path, help="Marker-panel YAML file.")
    parser.add_argument(
        "--he-source-root",
        type=Path,
        help="Parent of FD*/spatial.tar.gz; enables H&E overlays.",
    )
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def validate_filtered_input(adata, mouse, expected):
    """Require completed QC v1 raw counts without filtering a second time."""
    policy = adata.uns.get("qc_filter", {})
    if policy.get("version") != HD_QC_FILTER_VERSION or policy.get(
        "normalization_applied", True
    ):
        raise ValueError(f"{mouse}: expected unnormalized QC v1 input.")
    if (
        adata.uns.get("visium_hd", {}).get("bin_size_um") != 8
        or not adata.obs["in_tissue"].eq(1).all()
    ):
        raise ValueError(f"{mouse}: expected in-tissue 8 µm bins.")
    if adata.n_obs != int(expected["bins_retained"]) or adata.n_vars != int(
        expected["features_retained"]
    ):
        raise ValueError(f"{mouse}: matrix dimensions differ from filter_summary.csv.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError(f"{mouse}: bin and gene names must be unique.")
    if set(adata.obs["sample_id"]) != {mouse} or not adata.obs["qc_keep"].all():
        raise ValueError(f"{mouse}: sample IDs or retained-bin flags are inconsistent.")
    if not select_hd_qc_bins(adata.obs)["qc_keep"].all():
        raise ValueError(f"{mouse}: input contains bins failing the selected QC.")
    values = adata.X.data if hasattr(adata.X, "tocsr") else np.asarray(adata.X)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.equal(values, np.floor(values)).all()
    ):
        raise ValueError(
            f"{mouse}: expression must contain finite, nonnegative integer counts."
        )
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    if not np.array_equal(totals, adata.obs["total_counts"].to_numpy()):
        raise ValueError(f"{mouse}: matrix counts differ from input QC totals.")
    positions = adata.obsm["spatial"]
    if positions.shape != (adata.n_obs, 2) or not np.isfinite(positions).all():
        raise ValueError(f"{mouse}: invalid spatial coordinates.")


def cluster_hd_mouse(
    mouse, input_dir, output_dir, expected, markers, settings, he_root
):
    """Produce clustering and marker evidence; leave tissue labels for review."""
    source = input_dir / mouse / "filtered_counts.h5ad"
    adata = load_data(source)
    validate_filtered_input(adata, mouse, expected)
    image, scale = None, 1.0
    selection = None
    if he_root is not None:
        image, scale = load_visium_hd_histology(
            Path(adata.uns["visium_hd"]["source_path"]),
            he_root / mouse / "spatial.tar.gz",
        )
        selection = pd.read_parquet(input_dir / mouse / "bin_selection.parquet")
        if not selection.index.is_unique:
            raise ValueError(f"{mouse}: duplicate barcodes in bin_selection.parquet.")
        kept = selection.loc[selection["qc_keep"]]
        if len(kept) != adata.n_obs or not adata.obs_names.isin(kept.index).all():
            raise ValueError(f"{mouse}: applied mask does not match retained bins.")
        if not np.array_equal(
            kept.loc[adata.obs_names, ["x", "y"]], adata.obsm["spatial"]
        ):
            raise ValueError(f"{mouse}: applied mask and count coordinates differ.")

    sample_dir = output_dir / mouse
    sample_dir.mkdir()
    normalize_and_log(adata, target_sum=10000)
    cluster_and_umap(
        adata,
        **settings,
        exclude_genes=list(adata.uns["qc_filter"]["mt_genes_matched"]),
    )
    ranks = rank_annotation_markers(adata)
    ranks.to_csv(sample_dir / "cluster_markers.csv", index=False)
    coverage, expression = summarize_annotation_markers(adata, markers["panels"])
    coverage.to_csv(sample_dir / "marker_coverage.csv", index=False)
    expression.to_csv(sample_dir / "marker_expression.csv", index=False)
    summary = adata.obs.groupby("leiden", observed=True).agg(
        bins=("total_counts", "size"),
        median_umis=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"),
        median_mt_pct=("pct_counts_mt", "median"),
    )
    summary["mouse_id"] = mouse
    summary["top_genes"] = [
        "; ".join(
            ranks.loc[ranks["group"].astype(str).eq(str(cluster)), "names"].head(10)
        )
        for cluster in summary.index
    ]
    summary["tissue_annotation"] = ""
    summary["annotation_status"] = "unreviewed"
    summary["notes"] = ""
    summary.index.name = "cluster"
    summary.to_csv(sample_dir / "cluster_annotation.csv")
    plot_hd_annotation(
        adata,
        mouse,
        markers["panels"],
        markers["spatial_markers"],
        sample_dir,
        image=image,
        scale=scale,
    )
    if selection is not None:
        plot_hd_qc_histology(selection, mouse, image, scale, sample_dir, applied=True)
    adata.uns["tissue_annotation"] = {
        "status": "clusters_ready_for_review",
        "source": str(source),
        "normalization": "log2(1 + counts per 10000)",
        "raw_counts_layer": "counts",
        "batch_correction": False,
        "cluster_ids_are_sample_specific": True,
        "mt_excluded_from_pca": True,
        "settings": settings,
        "marker_ranking": "t-test_overestim_var; descriptive cluster-vs-rest only",
    }
    temporary = sample_dir / "clustered_bins.partial.h5ad"
    adata.write_h5ad(temporary, compression="gzip")
    temporary.replace(sample_dir / "clustered_bins.h5ad")
    return {
        "mouse_id": mouse,
        "bins": adata.n_obs,
        "features": adata.n_vars,
        "clusters": len(summary),
        "hvg_genes": int(adata.var["highly_variable"].sum()),
        "marker_genes_present": int(coverage["in_matrix"].sum()),
        "annotation_status": "unreviewed",
        "he_overlays": image is not None,
    }


def run_visium_hd_annotation(args, config):
    """Run each mouse independently on a compute node; save compact review files."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Run experimental-data annotation in a SLURM compute job.")
    # Check the optional clustering dependency before opening large matrices.
    try:
        importlib.metadata.version("igraph")
    except importlib.metadata.PackageNotFoundError as exc:
        raise ImportError(
            "Install annotation dependencies: pip install -e '.[annotation]'"
        ) from exc
    mice = [f"FD{i}" for i in range(1, 13)] if args.all_mice else args.mice
    if (
        not mice
        or len(mice) != len(set(mice))
        or any(re.fullmatch(r"FD(?:[1-9]|1[0-2])", mouse) is None for mouse in mice)
    ):
        raise ValueError("Specify --all-mice or unique --mice FD1 … FD12.")
    repo = get_config_path(config, "repo_root", ".")
    paths = get_config_section(config, "paths")
    input_dir = (
        args.input_dir
        or get_section_path(config, paths, "preprocess_output", "results/preprocess")
        / "visium_hd_008um"
        / HD_QC_FILTER_VERSION
    ).resolve()
    output_dir = (
        args.output_dir
        or get_section_path(config, paths, "integration_output", "results/integration")
        / "visium_hd_008um"
    ).resolve()
    for protected in [repo, input_dir, get_config_path(config, "data_dir", "data")]:
        if (
            output_dir == protected
            or protected in output_dir.parents
            or output_dir in protected.parents
        ):
            raise ValueError(
                f"Keep annotation results separate from inputs and code: {protected}"
            )
    if output_dir.exists():
        raise FileExistsError(
            f"Preserve existing results; choose an unused --output-dir: {output_dir}"
        )
    prior_run = json.loads((input_dir / "run.json").read_text())
    if prior_run["status"] != "complete" or not set(mice).issubset(
        prior_run["completed_mice"]
    ):
        raise ValueError("The selected QC run is incomplete for these mice.")
    expected = pd.read_csv(input_dir / "filter_summary.csv").set_index("mouse_id")
    if not expected.index.is_unique or not set(mice).issubset(expected.index):
        raise ValueError("Expected one filter summary row per requested mouse.")
    marker_path = args.markers or repo / "configs/thymus_markers.yaml"
    markers = yaml.safe_load(marker_path.read_text())
    if not isinstance(markers.get("panels"), dict) or not markers["panels"]:
        raise ValueError("Marker configuration must contain a nonempty panels mapping.")
    for genes in [*markers["panels"].values(), markers["spatial_markers"]]:
        if (
            not isinstance(genes, list)
            or not genes
            or not all(isinstance(gene, str) for gene in genes)
            or len(genes) != len(set(genes))
        ):
            raise ValueError("Each marker panel must be a list of unique gene symbols.")
    he_root = args.he_source_root.resolve() if args.he_source_root else None
    for mouse in mice:
        required = [input_dir / mouse / "filtered_counts.h5ad"]
        if he_root is not None:
            required += [
                input_dir / mouse / "bin_selection.parquet",
                he_root / mouse / "spatial.tar.gz",
            ]
        for path in required:
            if not path.is_file():
                raise FileNotFoundError(path)
    settings = {
        "n_top_genes": args.n_top_genes,
        "n_pcs": args.n_pcs,
        "n_neighbors": args.n_neighbors,
        "resolution": args.resolution,
        "random_state": args.seed,
    }
    run = {
        "status": "running",
        "mice": mice,
        "completed_mice": [],
        "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "settings": settings,
        "normalization": "log2(1 + counts per 10000)",
        "batch_correction": False,
        "markers": markers,
        "annotation_status": "unreviewed",
        "he_source_root": str(he_root) if he_root else None,
        "input_run_sha256": hashlib.sha256(
            (input_dir / "run.json").read_bytes()
        ).hexdigest(),
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["scanpy", "anndata", "igraph", "numpy", "scipy"]
        },
    }
    output_dir.mkdir(parents=True)
    run_path = output_dir / "run.json"
    run_path.write_text(json.dumps(run, indent=2) + "\n")
    started = time.monotonic()
    summaries = []
    try:
        for mouse in mice:
            print(
                f"{mouse}: clustering filtered bins for tissue annotation", flush=True
            )
            summaries.append(
                cluster_hd_mouse(
                    mouse,
                    input_dir,
                    output_dir,
                    expected.loc[mouse],
                    markers,
                    settings,
                    he_root,
                )
            )
            pd.DataFrame(summaries).to_csv(
                output_dir / "annotation_summary.csv", index=False
            )
            run["completed_mice"].append(mouse)
            run_path.write_text(json.dumps(run, indent=2) + "\n")
            gc.collect()
        run["status"] = "complete"
    except Exception as exc:
        run.update(status="failed", error=str(exc))
        raise
    finally:
        run["elapsed_seconds"] = time.monotonic() - started
        run_path.write_text(json.dumps(run, indent=2) + "\n")
    print(f"Annotation evidence ready for review: {output_dir}", flush=True)
    return 0


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    if args.visium_hd:
        return run_visium_hd_annotation(args, config)

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
