"""Tests for data loading helpers."""

import argparse
import hashlib
import io
import json
import runpy
import shutil
import tarfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from spatial_transcriptomics import data as data_module
from spatial_transcriptomics.data import (
    apply_hd_qc_filter,
    load_data,
    load_reference_genes,
    load_visium_hd_bin,
    load_visium_hd_histology,
    read_de_csv,
    stage_visium_hd_qc,
)


def _hd_filter_fixture():
    """Small raw-count input with threshold boundaries and a high-content bin."""
    ad = pytest.importorskip("anndata")
    from scipy import sparse

    counts = np.zeros((6, 6002), dtype=np.float32)
    counts[0, :20] = [5] + [1] * 18 + [2]  # Exactly 25 UMIs, 20 genes, 20% MT.
    counts[1, :20] = [4] + [1] * 18 + [2]  # UMI-only failure.
    counts[2, 1:20] = [7] + [1] * 18  # Gene-only failure.
    counts[3, :20] = [6] + [1] * 18 + [2]  # MT-only failure.
    counts[5, :6001] = 1  # Above both legacy upper cutoffs.
    counts[5, 1] = 40000
    adata = ad.AnnData(
        sparse.csr_matrix(counts),
        obs=pd.DataFrame(
            {"in_tissue": [1] * 6, "pct_counts_rp": [99.0] * 6},
            index=[f"bin-{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["mt-Nd1"] + [f"Gene{i}" for i in range(6001)]),
    )
    adata.obsm["spatial"] = np.arange(12, dtype=float).reshape(6, 2)
    adata.uns["visium_hd"] = {"bin_size_um": 8}
    return adata


def test_apply_hd_filter_preserves_counts_features_coordinates_and_input(tmp_path):
    pytest.importorskip("scanpy")
    import anndata as ad

    source = _hd_filter_fixture()
    original = source.copy()
    kept, selection = apply_hd_qc_filter(source, ["mt-Nd1"])
    assert selection["qc_keep"].tolist() == [True, False, False, False, False, True]
    assert kept.obs_names.tolist() == ["bin-0", "bin-5"]
    assert kept.n_vars == 6002  # Includes an all-zero gene.
    assert (kept.X != original.X[[0, 5]]).nnz == 0
    np.testing.assert_array_equal(
        kept.obsm["spatial"], original.obsm["spatial"][[0, 5]]
    )
    pd.testing.assert_frame_equal(source.obs, original.obs)
    assert (source.X != original.X).nnz == 0
    assert "qc_filter" not in source.uns
    kept.write_h5ad(tmp_path / "filtered.h5ad")
    restored = ad.read_h5ad(tmp_path / "filtered.h5ad")
    assert (restored.X != kept.X).nnz == 0
    assert not restored.uns["qc_filter"]["normalization_applied"]
    with pytest.raises(ValueError, match="mitochondrial"):
        apply_hd_qc_filter(source, ["not_present"])
    source.obs.loc["bin-0", "in_tissue"] = 0
    with pytest.raises(ValueError, match="in-tissue"):
        apply_hd_qc_filter(source, ["mt-Nd1"])


def test_hd_filter_job_roundtrip_and_review_consistency(tmp_path, monkeypatch):
    """Real sparse H5 → archive → filter → H5AD, with unchanged source inputs."""
    h5py = pytest.importorskip("h5py")
    import anndata as ad

    repo_root = Path(__file__).resolve().parents[1]
    script = runpy.run_path(str(repo_root / "scripts" / "02_preprocess.py"))
    source = _hd_filter_fixture()
    fixture = tmp_path / "fixture" / "square_008um"
    (fixture / "spatial").mkdir(parents=True)
    matrix = source.X.T.tocsc()
    with h5py.File(fixture / "filtered_feature_bc_matrix.h5", "w") as handle:
        group = handle.create_group("matrix")
        for name, values in {
            "data": matrix.data,
            "indices": matrix.indices,
            "indptr": matrix.indptr,
            "shape": matrix.shape,
            "barcodes": np.array(source.obs_names, dtype="S"),
        }.items():
            group.create_dataset(name, data=values)
        features = group.create_group("features")
        for name, values in {
            "id": [f"ENSM{i}" for i in range(source.n_vars)],
            "name": source.var_names,
            "feature_type": ["Gene Expression"] * source.n_vars,
            "genome": ["mm10"] * source.n_vars,
        }.items():
            features.create_dataset(name, data=np.array(values, dtype="S"))
    pd.DataFrame(
        {
            "barcode": source.obs_names,
            "in_tissue": 1,
            "pxl_col_in_fullres": source.obsm["spatial"][:, 0],
            "pxl_row_in_fullres": source.obsm["spatial"][:, 1],
        }
    ).iloc[::-1].to_parquet(fixture / "spatial" / "tissue_positions.parquet")
    (fixture / "spatial" / "scalefactors_json.json").write_text(
        '{"tissue_hires_scalef": 0.5}'
    )
    _, selection = apply_hd_qc_filter(source, ["mt-Nd1"])
    for mouse in ["FD1", "FD2"]:
        archive_path = tmp_path / "raw" / mouse / "binned_outputs.tar.gz"
        archive_path.parent.mkdir(parents=True)
        with tarfile.open(archive_path, "w:gz") as archive:
            archive.add(fixture, arcname="binned_outputs/square_008um")
        qc_path = tmp_path / "qc" / mouse / "bin_qc.parquet"
        qc_path.parent.mkdir(parents=True)
        saved = selection.iloc[::-1].copy()  # Barcode alignment, not row position.
        saved.loc["bin-0", "pct_counts_mt"] = 20.000000298023224  # Old float32 result.
        saved.to_parquet(qc_path)
    originals = {
        p: hashlib.sha256(p.read_bytes()).hexdigest()
        for root in [tmp_path / "raw", tmp_path / "qc"]
        for p in root.rglob("*")
        if p.is_file()
    }
    args = argparse.Namespace(
        config="synthetic fixture",
        all_mice=False,
        mice=["FD1", "FD2"],
        source_root=tmp_path / "raw",
        qc_dir=tmp_path / "qc",
        output_dir=tmp_path / "filtered",
        scratch_dir=tmp_path / "scratch",
    )
    config = {"repo_root": str(repo_root), "visium_hd": {"bin_size_um": 8}}
    run_filter = script["run_visium_hd_filter"]
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="SLURM"):
        run_filter(config, args)
    monkeypatch.setenv("SLURM_JOB_ID", "synthetic-fixture-only")
    assert run_filter(config, args) == 0
    run = json.loads((args.output_dir / "run.json").read_text())
    assert run["status"] == "complete"
    assert run["completed_mice"] == ["FD1", "FD2"]
    result = ad.read_h5ad(args.output_dir / "FD1" / "filtered_counts.h5ad")
    assert result.obs_names.tolist() == ["bin-0", "bin-5"]
    assert (result.X != source.X[[0, 5]]).nnz == 0
    assert result.uns["spatial"]["FD1"]["scalefactors"]["tissue_hires_scalef"] == 0.5
    assert result.obs["sample_id"].eq("FD1").all()
    decisions = pd.read_parquet(args.output_dir / "FD1" / "bin_selection.parquet")
    assert len(decisions) == source.n_obs
    summary = pd.read_csv(args.output_dir / "filter_summary.csv")
    assert summary["bins_retained"].tolist() == [2, 2]
    assert summary["bins_excluded"].tolist() == [4, 4]
    assert summary["mt_boundary_corrections"].tolist() == [1, 1]
    assert len(pd.read_csv(args.output_dir / "FD1" / "failure_summary.csv")) == 11
    assert not list(args.scratch_dir.iterdir())
    for path, digest in originals.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    with pytest.raises(FileExistsError, match="Preserve existing"):
        run_filter(config, args)

    # A changed matrix/review cannot silently create a second 'successful' result.
    args.output_dir = tmp_path / "mismatch"
    changed = selection.copy()
    changed.loc["bin-0", "total_counts"] = 999
    changed.to_parquet(args.qc_dir / "FD1" / "bin_qc.parquet")
    with pytest.raises(ValueError, match="differs from the reviewed QC"):
        run_filter(config, args)
    assert json.loads((args.output_dir / "run.json").read_text())["status"] == "failed"
    assert not list(args.output_dir.rglob("filtered_counts.h5ad"))


def test_load_reference_genes_returns_nonempty_lines(tmp_path: Path) -> None:
    """Gene lists should ignore blank lines."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("GeneA\n\nGeneB\n")

    assert load_reference_genes(gene_file) == ["GeneA", "GeneB"]


def test_load_reference_genes_missing_file_returns_empty_list(tmp_path: Path) -> None:
    """Missing gene files should produce an empty list."""
    assert load_reference_genes(tmp_path / "missing.txt") == []


def test_load_data_raises_for_missing_path(tmp_path: Path) -> None:
    """Missing data paths should raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_data(tmp_path / "missing.h5ad")


def test_load_visium_hd_bin_attaches_coordinates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The HD loader should join parquet positions in matrix-barcode order."""
    bin_path = tmp_path / "binned_outputs" / "square_008um"
    spatial_path = bin_path / "spatial"
    spatial_path.mkdir(parents=True)
    (bin_path / "filtered_feature_bc_matrix.h5").touch()
    (spatial_path / "tissue_positions.parquet").touch()

    adata = SimpleNamespace(
        obs=pd.DataFrame(index=pd.Index(["bin-b", "bin-a"])),
        obs_names=pd.Index(["bin-b", "bin-a"]),
        var_names=pd.Index(["GeneA"]),
        obsm={},
        uns={},
        var_names_make_unique=lambda: None,
    )
    fake_scanpy = SimpleNamespace(read_10x_h5=lambda *_args, **_kwargs: adata)
    monkeypatch.setattr(data_module, "_load_scanpy", lambda: fake_scanpy)
    monkeypatch.setattr(
        data_module.pd,
        "read_parquet",
        lambda _path: pd.DataFrame(
            {
                "barcode": ["bin-a", "bin-b"],
                "in_tissue": [1, 1],
                "array_row": [10, 20],
                "array_col": [30, 40],
                "pxl_row_in_fullres": [100.0, 200.0],
                "pxl_col_in_fullres": [300.0, 400.0],
            },
        ),
    )

    result = load_visium_hd_bin(
        tmp_path,
        bin_size_um=8,
        library_id="FD_1",
        load_images=False,
    )

    np.testing.assert_array_equal(
        result.obsm["spatial"],
        np.array([[400.0, 200.0], [300.0, 100.0]]),
    )
    assert result.obs["bin_size_um"].tolist() == [8, 8]
    assert result.uns["visium_hd"]["library_id"] == "FD_1"


def test_read_de_csv_reads_valid_file(tmp_path: Path) -> None:
    """Valid DE CSV files should be loaded unchanged."""
    csv_file = tmp_path / "de.csv"
    pd.DataFrame(
        {"stat": [1.2, 0.5], "padj": [0.01, 0.2]},
        index=["GeneA", "GeneB"],
    ).to_csv(csv_file)

    result = read_de_csv(csv_file)

    assert list(result.index) == ["GeneA", "GeneB"]
    assert list(result.columns) == ["stat", "padj"]


def test_read_de_csv_raises_for_missing_required_columns(tmp_path: Path) -> None:
    """Missing required DE columns should raise a clear error."""
    csv_file = tmp_path / "bad.csv"
    pd.DataFrame({"stat": [1.2]}, index=["GeneA"]).to_csv(csv_file)

    with pytest.raises(ValueError, match="missing columns"):
        read_de_csv(csv_file)


def test_read_de_csv_drops_invalid_rows_and_duplicate_genes(tmp_path: Path) -> None:
    """Infinite values, missing values, and duplicate genes should be removed."""
    csv_file = tmp_path / "filtered.csv"
    pd.DataFrame(
        {
            "stat": [1.2, float("inf"), 0.7, 0.9],
            "padj": [0.01, 0.02, None, 0.03],
        },
        index=["GeneA", "GeneB", "GeneC", "GeneA"],
    ).to_csv(csv_file)

    result = read_de_csv(csv_file)

    assert list(result.index) == ["GeneA"]
    assert result.loc["GeneA", "stat"] == 1.2


def test_hd_staging_extracts_only_requested_inputs(tmp_path: Path) -> None:
    """Archive prefixes cannot escape scratch; other resolutions/images stay out."""
    archive = tmp_path / "input.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        for name in [
            "prefix/square_002um/filtered_feature_bc_matrix.h5",
            "../../square_008um/filtered_feature_bc_matrix.h5",
            "prefix/square_008um/spatial/tissue_hires_image.png",
            "prefix/square_008um/spatial/tissue_positions.parquet",
        ]:
            info = tarfile.TarInfo(name)
            info.size = 3
            handle.addfile(info, io.BytesIO(b"abc"))
    staged = stage_visium_hd_qc(archive, tmp_path / "scratch")
    assert (staged / "filtered_feature_bc_matrix.h5").read_bytes() == b"abc"
    assert (staged / "spatial/tissue_positions.parquet").read_bytes() == b"abc"
    assert (
        len([path for path in (tmp_path / "scratch").rglob("*") if path.is_file()]) == 2
    )
    with pytest.raises(FileNotFoundError, match="missing"):
        stage_visium_hd_qc(archive, tmp_path / "missing", bin_size_um=16)


def test_hd_histology_reads_image_and_matching_bin_scale(tmp_path: Path) -> None:
    """Read separate archives, select 8 µm metadata, and preserve image pixels."""
    import matplotlib.image as mpimg

    original = np.zeros((12, 16, 3), dtype=np.uint8)
    original[2:5, 7:9] = [255, 0, 128]
    image_buffer = io.BytesIO()
    mpimg.imsave(image_buffer, original, format="png")
    binned = tmp_path / "binned_outputs.tar.gz"
    spatial = tmp_path / "spatial.tar.gz"
    for archive_path, entries in [
        (
            binned,
            {
                "square_002um/spatial/scalefactors_json.json": json.dumps(
                    {"tissue_hires_scalef": 1}
                ).encode(),
                "prefix/square_008um/spatial/scalefactors_json.json": json.dumps(
                    {"tissue_hires_scalef": 0.25}
                ).encode(),
            },
        ),
        (spatial, {"tissue_hires_image.png": image_buffer.getvalue()}),
    ]:
        with tarfile.open(archive_path, "w:gz") as archive:
            for name, content in entries.items():
                info = tarfile.TarInfo(name)
                info.size = len(content)
                archive.addfile(info, io.BytesIO(content))
    image, scale = load_visium_hd_histology(binned, spatial)
    assert scale == 0.25
    np.testing.assert_allclose(image[:, :, :3], original / 255, atol=1e-7)
    with pytest.raises(FileNotFoundError, match="tissue_hires_image"):
        load_visium_hd_histology(binned, binned)


def test_hd_batch_outputs_and_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise sparse loading, QC, plots, cache relocation, and invalidation."""
    h5py = pytest.importorskip("h5py")
    pytest.importorskip("scanpy")
    pytest.importorskip("pyarrow")
    from scipy import sparse

    repo_root = Path(__file__).resolve().parents[1]
    script = runpy.run_path(str(repo_root / "scripts" / "01_qc.py"))
    run_qc = script["run_visium_hd_qc"]
    refs = tmp_path / "reference_genomes"
    refs.mkdir()
    (refs / "mouse_mitochondrial_genes.txt").write_text("mt-Test\n")
    (refs / "mouse_ribosomal_genes.txt").write_text("RplTest\n")
    counts = np.array([[0, 0, 0], [2, 3, 5], [5, 10, 10], [10, 20, 20]])
    barcodes = [
        "s_008um_00000_00000-1",
        "s_008um_00000_00001-1",
        "s_008um_00001_00000-1",
        "s_008um_00001_00001-1",
    ]
    fixture = tmp_path / "fixture" / "binned_outputs" / "square_008um"
    (fixture / "spatial").mkdir(parents=True)
    (fixture / "spatial" / "scalefactors_json.json").write_text(
        '{"tissue_hires_scalef": 0.5}'
    )
    matrix = sparse.csc_matrix(counts.T)
    with h5py.File(fixture / "filtered_feature_bc_matrix.h5", "w") as handle:
        group = handle.create_group("matrix")
        for name, values in {
            "data": matrix.data,
            "indices": matrix.indices,
            "indptr": matrix.indptr,
            "shape": matrix.shape,
            "barcodes": np.array(barcodes, dtype="S"),
        }.items():
            group.create_dataset(name, data=values)
        features = group.create_group("features")
        for name, values in {
            "id": ["ENSM1", "ENSM2", "ENSM3"],
            "name": ["mt-Test", "RplTest", "GeneA"],
            "feature_type": ["Gene Expression"] * 3,
            "genome": ["mm10"] * 3,
        }.items():
            features.create_dataset(name, data=np.array(values, dtype="S"))
    pd.DataFrame(
        {
            "barcode": barcodes[::-1],
            "pxl_col_in_fullres": [20, 10, 20, 10],
            "pxl_row_in_fullres": [20, 20, 10, 10],
            "in_tissue": [1] * 4,
        }
    ).to_parquet(fixture / "spatial" / "tissue_positions.parquet")
    source = tmp_path / "source"
    for mouse in ["FD1", "FD2"]:
        (source / mouse).mkdir(parents=True)
        with tarfile.open(source / mouse / "binned_outputs.tar.gz", "w:gz") as handle:
            handle.add(fixture, arcname="binned_outputs/square_008um")
    output = tmp_path / "output"
    scratch = tmp_path / "scratch"
    args = argparse.Namespace(
        config="synthetic",
        mice=["FD1", "FD2"],
        all_mice=False,
        source_root=source,
        scratch_dir=scratch,
        output_dir=output,
        force=False,
        no_plots=False,
        he_overlays=False,
        he_image_root=None,
    )
    config = {"repo_root": str(tmp_path), "visium_hd": {"bin_size_um": 8}}
    assert run_qc(config, args) == 0
    qc = pd.read_parquet(output / "FD1" / "bin_qc.parquet")
    np.testing.assert_allclose(qc["total_counts"], [0, 10, 25, 50])
    np.testing.assert_allclose(qc["pct_counts_mt"], [0, 20, 20, 20])
    np.testing.assert_allclose(qc["pct_counts_rp"], [0, 30, 40, 40])
    assert qc["x"].tolist() == [10, 20, 10, 20]
    assert qc["y"].tolist() == [10, 10, 20, 20]
    comparison = output / "comparisons" / "FD1_FD2"
    summary = pd.read_csv(comparison / "cross_mouse_qc_summary.csv")
    assert summary["bins"].tolist() == [4, 4]
    failures = pd.read_csv(comparison / "cross_mouse_failure_summary.csv")
    exclusive = failures[failures["group_type"] == "exclusive"]
    assert exclusive.groupby(["mouse_id", "candidate"])["bins"].sum().eq(4).all()
    assert exclusive.groupby(["mouse_id", "candidate"])["umis"].sum().eq(85).all()
    pd.testing.assert_frame_equal(
        pd.read_csv(output / "FD1" / "failure_summary.csv"),
        failures[failures["mouse_id"] == "FD1"].reset_index(drop=True),
    )
    assert (comparison / "cross_mouse_qc_distributions.png").is_file()
    assert len(list(output.rglob("*.png"))) == 9
    assert not list(scratch.iterdir())

    def fail_load(*_args, **_kwargs):
        raise AssertionError("Unexpected matrix reload")

    monkeypatch.setitem(run_qc.__globals__, "load_visium_hd_bin", fail_load)
    args.no_plots = True
    (output / "FD1" / "failure_summary.csv").unlink()
    assert run_qc(config, args) == 0
    assert (output / "FD1" / "failure_summary.csv").exists()
    pd.testing.assert_frame_equal(
        pd.read_csv(comparison / "cross_mouse_failure_summary.csv"),
        failures,
    )

    # Copying results to a new configured root must reuse the cached metrics.
    previous_output = output
    output = tmp_path / "new_results" / "qc" / "visium_hd_008um"
    shutil.copytree(previous_output, output)
    cached_table = (output / "FD1" / "bin_qc.parquet").read_bytes()
    config["paths"] = {"qc_output": str(output.parent)}
    config["visium_hd"] = {"bin_size_um": 8, "source_dir": str(source / "FD1")}
    args.source_root = None
    args.output_dir = None
    assert run_qc(config, args) == 0
    assert (output / "FD1" / "bin_qc.parquet").read_bytes() == cached_table
    pd.testing.assert_frame_equal(
        pd.read_csv(output / "comparisons" / "FD1_FD2" / "cross_mouse_qc_summary.csv"),
        summary,
    )
    assert (previous_output / "FD1" / "bin_qc.parquet").read_bytes() == cached_table

    # H&E can live separately; overlay generation must not reload expression data.
    import matplotlib.image as mpimg

    image_buffer = io.BytesIO()
    mpimg.imsave(image_buffer, np.ones((20, 20, 3)), format="png")
    image_root = tmp_path / "images"
    for mouse in args.mice:
        (image_root / mouse).mkdir(parents=True)
        with tarfile.open(image_root / mouse / "spatial.tar.gz", "w:gz") as archive:
            content = image_buffer.getvalue()
            info = tarfile.TarInfo("tissue_hires_image.png")
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    args.he_image_root = image_root
    args.he_overlays = True
    with pytest.raises(ValueError, match="--no-plots"):
        run_qc(config, args)
    args.no_plots = False
    assert run_qc(config, args) == 0
    for mouse in args.mice:
        assert (output / mouse / "qc_failures_he.png").is_file()
    assert (output / "FD1" / "bin_qc.parquet").read_bytes() == cached_table
    args.he_overlays = False
    args.no_plots = True

    # An interrupted marker write must trigger recomputation, including --force.
    metadata_path = output / "FD1" / "qc_metadata.json"
    saved_metadata = metadata_path.read_text()
    metadata_path.write_text('{"input":')
    for force in [False, True]:
        args.force = force
        with pytest.raises(AssertionError, match="Unexpected matrix reload"):
            run_qc(config, args)
    args.force = False
    metadata_path.write_text(saved_metadata)

    # Changing a reference list must invalidate the cache.
    (refs / "mouse_mitochondrial_genes.txt").write_text("mt-Test\nmt-New\n")
    with pytest.raises(AssertionError, match="Unexpected matrix reload"):
        run_qc(config, args)
