"""Tests for plotting helpers."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import QuadMesh

from spatial_transcriptomics.analysis import evaluate_hd_qc
from spatial_transcriptomics.plotting import (
    plot_hd_qc,
    plot_hd_qc_comparison,
    plot_hd_qc_histology,
    plot_hd_spatial,
    pretty_title,
    stacked_bar,
    summary_dotplot,
)

mpl.use("Agg")


def _plotting_frame() -> pd.DataFrame:
    """Build a small categorical dataframe for plotting tests."""
    frame = pd.DataFrame(
        {
            "sample": ["s1", "s1", "s2", "s2"],
            "condition": ["flash", "flash", "conv", "conv"],
            "celltype": ["T", "B", "T", "B"],
            "fraction": [0.6, 0.4, 0.7, 0.3],
            "celltype_is_immune": [True, True, True, True],
        },
    )
    frame["condition"] = pd.Categorical(
        frame["condition"],
        categories=["flash", "conv"],
    )
    return frame


def test_stacked_bar_writes_output_file(tmp_path: Path) -> None:
    """Stacked bar plots should be saved to disk."""
    outfile = tmp_path / "stacked.png"

    stacked_bar(_plotting_frame(), "fraction", "Fractions", outfile)

    assert outfile.exists()
    assert outfile.stat().st_size > 0


def test_summary_dotplot_writes_output_file(tmp_path: Path) -> None:
    """Summary dotplots should be saved to disk."""
    outfile = tmp_path / "dotplot.png"

    summary_dotplot(_plotting_frame(), "fraction", "Summary", outfile)

    assert outfile.exists()
    assert outfile.stat().st_size > 0


def test_pretty_title_formats_known_conditions() -> None:
    """Known contrast labels should be expanded for presentation."""
    assert pretty_title("flash_vs_conv") == "FLASH vs Conventional"


def _hd_plot_frame() -> pd.DataFrame:
    """Build shuffled, rotated bins with one absent position inside the grid."""
    rows = np.array([4, 3, 4, 3, 4])
    cols = np.array([8, 7, 9, 9, 7])
    return pd.DataFrame(
        {
            "x": 100 + 4 * cols - 3 * rows,
            "y": 200 + 3 * cols + 4 * rows,
            "total_counts": [80, 10, 90, 30, 50],
        },
        index=[f"s_008um_{row:05d}_{col:05d}-1" for row, col in zip(rows, cols)],
    )


@pytest.mark.parametrize("all_pass", [False, True])
def test_hd_histology_alignment_and_exclusive_failures(tmp_path, monkeypatch, all_pass):
    """Keep H&E orientation, transparent background, and threshold equality correct."""
    qc = _hd_plot_frame()
    qc["total_counts"] = [80, 10, 90, 25, 50]
    qc["n_genes_by_counts"] = [50, 8, 50, 20, 10]
    qc["pct_counts_mt"] = [21, 20, 20, 20, 21]
    if all_pass:
        qc[["total_counts", "n_genes_by_counts", "pct_counts_mt"]] = [25, 20, 20]
    original = qc.copy(deep=True)
    image = np.ones((150, 100, 3))
    inspected = []

    def inspect_figure(fig, outfile):
        assert outfile.name == "qc_failures_he.png"
        for ax in fig.axes:
            np.testing.assert_array_equal(ax.images[0].get_array(), image)
            assert ax.get_xlim() == (-0.5, 99.5)
            assert ax.get_ylim() == (149.5, -0.5)
        for ax, selected in zip(fig.axes[1:], [0, 1, 4]):
            background, overlay = ax.collections
            assert background.get_array().count() == 0
            assert overlay.get_array().count() == (0 if all_pass else 1)
            if not all_pass:
                corners = overlay.get_coordinates()
                centers = (corners[:-1, :-1] + corners[1:, 1:]) / 2
                mask = ~np.ma.getmaskarray(overlay.get_array()).reshape(
                    centers.shape[:2]
                )
                np.testing.assert_allclose(
                    centers[mask][0], qc.iloc[selected][["x", "y"]].to_numpy() * 0.5
                )
            assert overlay.get_alpha() == 0.55
        fig.savefig(outfile)
        inspected.append(True)

    monkeypatch.setattr(
        "spatial_transcriptomics.plotting._save_qc_figure", inspect_figure
    )
    plot_hd_qc_histology(qc, "synthetic", image, 0.5, tmp_path)
    assert inspected == [True]
    pd.testing.assert_frame_equal(qc, original)


def test_hd_spatial_preserves_values_masks_and_affine_geometry() -> None:
    """Filled bins retain actual values, orientation, gaps, and rejected bins."""
    qc = _hd_plot_frame()
    mask = np.array([True, True, False, True, False])
    values = qc["total_counts"].to_numpy(copy=True)
    original_qc = qc.copy(deep=True)
    original_mask = mask.copy()
    original_values = values.copy()
    fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
    try:
        mesh = plot_hd_spatial(ax, qc, mask, values)
        assert isinstance(mesh, QuadMesh)
        rendered_values = mesh.get_array()
        assert isinstance(rendered_values, np.ma.MaskedArray)
        np.testing.assert_array_equal(
            np.sort(rendered_values.compressed()), [10, 30, 80]
        )
        assert mesh.get_clim() == (10, 90)
        assert ax.yaxis_inverted()
        assert ax.get_aspect() == 1

        # Grid corners must follow the rotated array, not an axis-aligned scatter.
        cols, rows = np.meshgrid(np.arange(6.5, 10), np.arange(2.5, 5))
        expected_corners = np.stack(
            [100 + 4 * cols - 3 * rows, 200 + 3 * cols + 4 * rows], axis=-1
        )
        np.testing.assert_allclose(
            np.asarray(mesh.get_coordinates(), dtype=float), expected_corners
        )

        canvas = fig.canvas
        assert isinstance(canvas, FigureCanvasAgg)
        canvas.draw()
        image = np.asarray(canvas.buffer_rgba())

        def rendered_color(row: int, col: int) -> np.ndarray:
            point = (100 + 4 * col - 3 * row, 200 + 3 * col + 4 * row)
            x, y = ax.transData.transform(point)
            return image[image.shape[0] - 1 - int(y), int(x), :3]

        np.testing.assert_array_equal(rendered_color(3, 8), [255, 255, 255])
        np.testing.assert_allclose(rendered_color(4, 9), [211, 211, 211], atol=1)
        assert mesh.norm is not None
        assert mesh.cmap is not None
        expected_color = np.array(mesh.cmap(mesh.norm(80)))[:3] * 255
        np.testing.assert_allclose(rendered_color(4, 8), expected_color, atol=1)
    finally:
        plt.close(fig)

    pd.testing.assert_frame_equal(qc, original_qc)
    np.testing.assert_array_equal(mask, original_mask)
    np.testing.assert_array_equal(values, original_values)


def test_hd_spatial_empty_selection_keeps_shared_color_limits() -> None:
    """A threshold retaining no bins must still use the full sample's scale."""
    qc = _hd_plot_frame()
    fig, ax = plt.subplots()
    try:
        mesh = plot_hd_spatial(
            ax, qc, np.zeros(len(qc), dtype=bool), qc["total_counts"].to_numpy()
        )
        rendered_values = mesh.get_array()
        assert isinstance(rendered_values, np.ma.MaskedArray)
        assert rendered_values.count() == 0
        assert mesh.get_clim() == (10, 90)
    finally:
        plt.close(fig)


def test_hd_qc_includes_unfiltered_mt_on_its_percentage_scale(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The MT spatial panel must retain percentages, including high-MT bins."""
    qc = _hd_plot_frame().assign(
        n_genes_by_counts=[70, 9, 80, 25, 40],
        pct_counts_mt=[0, 10, 20, 50, 100],
    )
    original = qc.copy(deep=True)
    saved = []

    def inspect_figure(fig, path):
        saved.append(path.name)
        if path.name == "qc_distributions.png":
            assert len(fig.axes) == 3
            assert fig.axes[2].get_xlabel() == "pct_counts_mt"
        if path.name == "spatial_qc.png":
            ax = next(
                ax for ax in fig.axes if ax.get_title() == "Mitochondrial counts (%)"
            )
            mesh = ax.collections[-1]
            assert mesh.get_clim() == (0, 100)
            np.testing.assert_array_equal(
                np.sort(mesh.get_array().compressed()),
                [0, 10, 20, 50, 100],
            )
        plt.close(fig)

    monkeypatch.setattr(
        "spatial_transcriptomics.plotting._save_qc_figure",
        inspect_figure,
    )
    plot_hd_qc(qc, "FD1", tmp_path)
    assert len(saved) == 3
    pd.testing.assert_frame_equal(qc, original)


@pytest.mark.parametrize("n_mice", [1, 12])
def test_combined_qc_distributions_preserve_all_bins(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, n_mice: int
) -> None:
    """Every mouse keeps its unfiltered histograms and original metric scales."""
    frames = {}
    qc_paths = {}
    candidates = []
    for i in range(1, n_mice + 1):
        mouse = f"FD{i}"
        qc = _hd_plot_frame().assign(
            n_genes_by_counts=[70, 9, 80, 25, 40],
            pct_counts_mt=[0, 10, 20, 50, 100],
        ).iloc[: 4 + i % 2].copy()
        qc["total_counts"] *= i
        frames[mouse] = qc
        qc_paths[mouse] = tmp_path / f"{mouse}.parquet"
        qc.to_parquet(qc_paths[mouse])
        summary, _ = evaluate_hd_qc(qc)
        candidates.append(summary.assign(mouse_id=mouse))
    saved = []

    def inspect_figure(fig, path):
        saved.append(path.name)
        if path.name == "cross_mouse_qc_distributions.png":
            assert len(fig.axes) == 3 * n_mice
            for i, (mouse, qc) in enumerate(frames.items()):
                row = fig.axes[3 * i : 3 * i + 3]
                assert row[0].get_ylabel().startswith(f"{mouse}\n")
                for ax, column in zip(
                    row, ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
                ):
                    values = qc[column].to_numpy()
                    if column != "pct_counts_mt":
                        values = np.log10(values + 1)
                    heights, edges = np.histogram(values, bins=80)
                    np.testing.assert_array_equal(
                        [patch.get_height() for patch in ax.patches], heights
                    )
                    np.testing.assert_allclose(
                        [patch.get_x() for patch in ax.patches], edges[:-1]
                    )
                    assert sum(heights) == len(qc)
        plt.close(fig)

    monkeypatch.setattr(
        "spatial_transcriptomics.plotting._save_qc_figure", inspect_figure
    )
    plot_hd_qc_comparison(pd.concat(candidates), qc_paths, tmp_path)
    assert saved.count("cross_mouse_qc_distributions.png") == 1


@pytest.mark.parametrize(
    "invalid_layout", ["barcode", "duplicate", "collinear", "warped"]
)
def test_hd_spatial_rejects_invalid_layout(invalid_layout: str) -> None:
    """Invalid or non-affine coordinates must not silently relocate expression."""
    qc = _hd_plot_frame()
    if invalid_layout == "barcode":
        qc.index = pd.Index(["unrecognized"] + qc.index.tolist()[1:])
    elif invalid_layout == "duplicate":
        qc.index = pd.Index([qc.index[1].replace("-1", "-2")] + qc.index.tolist()[1:])
    elif invalid_layout == "collinear":
        qc["y"] = qc["x"] * 2
    else:
        qc.iloc[0, 0] = qc["x"].to_numpy()[0] + 10

    fig, ax = plt.subplots()
    try:
        with pytest.raises(ValueError):
            plot_hd_spatial(ax, qc, np.ones(len(qc), dtype=bool))
    finally:
        plt.close(fig)
