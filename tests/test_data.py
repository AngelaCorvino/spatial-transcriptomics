"""Tests for data loading helpers."""

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
