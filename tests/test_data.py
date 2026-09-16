"""Tests for data loading helpers."""

import argparse
import io
import runpy
import tarfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from spatial_transcriptomics import data as data_module
from spatial_transcriptomics.data import (
    load_data,
    load_reference_genes,
    load_visium_hd_bin,
    read_de_csv,
    stage_visium_hd_qc,
)


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


def test_hd_batch_outputs_and_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise real sparse 10x loading, QC, plots, cache reuse, and invalidation."""
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
    barcodes = ["a", "b", "c", "d"]
    fixture = tmp_path / "fixture" / "binned_outputs" / "square_008um"
    (fixture / "spatial").mkdir(parents=True)
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
            "pxl_col_in_fullres": [40, 30, 20, 10],
            "pxl_row_in_fullres": [4, 3, 2, 1],
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
    )
    config = {"repo_root": str(tmp_path), "visium_hd": {"bin_size_um": 8}}
    assert run_qc(config, args) == 0
    qc = pd.read_parquet(output / "FD1" / "bin_qc.parquet")
    np.testing.assert_allclose(qc["total_counts"], [0, 10, 25, 50])
    np.testing.assert_allclose(qc["pct_counts_mt"], [0, 20, 20, 20])
    np.testing.assert_allclose(qc["pct_counts_rp"], [0, 30, 40, 40])
    assert qc["x"].tolist() == [10, 20, 30, 40]
    comparison = output / "comparisons" / "FD1_FD2"
    summary = pd.read_csv(comparison / "cross_mouse_qc_summary.csv")
    assert summary["bins"].tolist() == [4, 4]
    assert len(list(output.rglob("*.png"))) == 8
    assert not list(scratch.iterdir())

    def fail_load(*_args, **_kwargs):
        raise AssertionError("Unexpected matrix reload")

    monkeypatch.setitem(run_qc.__globals__, "load_visium_hd_bin", fail_load)
    args.no_plots = True
    assert run_qc(config, args) == 0

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
