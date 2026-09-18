"""Plotting utilities for spatial transcriptomics analysis.

Provides functions for visualization of cell type composition,
differential expression results, and other analysis outputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from spatial_transcriptomics.analysis import HD_QC_CANDIDATES, evaluate_hd_qc


def _save_qc_figure(fig, outfile: Path) -> None:
    """Save and immediately release a batch QC figure."""
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_hd_spatial(ax, qc: pd.DataFrame, mask, values=None, *, show_excluded=True):
    """Draw filled HD bin footprints without scatter-marker moire.

    ``qc`` contains image coordinates ``x``/``y`` and HD barcodes as its index.
    Barcode rows/columns define adjacent square bins; an affine fit preserves
    their image orientation. Reject layouts deviating by more than 0.1 bin.
    Missing bins stay transparent; excluded bins are gray unless hidden. Values and
    masks follow dataframe order and are never modified, smoothed or aggregated.
    Color limits use all supplied values, keeping candidate panels comparable.
    """
    positions = qc.index.to_series().str.extract(r"^s_\d+um_(\d+)_(\d+)-\d+$")
    if qc.empty or positions.isna().any().any():
        raise ValueError("Spatial HD plots require row/column Visium HD barcodes.")
    if positions.duplicated().any():
        raise ValueError("Spatial HD plots require unique bin positions.")
    row, col = positions.to_numpy(dtype=int).T
    row, col = row - row.min(), col - col.min()
    design = np.column_stack([col, row, np.ones(len(qc))])
    centers = qc[["x", "y"]].to_numpy(dtype=float)
    if not np.isfinite(centers).all():
        raise ValueError("Spatial HD plots require finite image coordinates.")
    transform, _, rank, _ = np.linalg.lstsq(design, centers, rcond=None)
    pitch = np.linalg.norm(transform[:2], axis=1).min()
    fitted = np.einsum("ij,jk->ik", design, transform)
    error = np.linalg.norm(fitted - centers, axis=1).max()
    if (
        rank < 3
        or pitch <= 0
        or abs(np.linalg.det(transform[:2])) < 1e-6 * pitch**2
        or error > 0.1 * pitch
    ):
        raise ValueError("Image coordinates do not match an affine HD bin grid.")

    mask = np.asarray(mask, dtype=bool)
    colors = np.ones(len(qc)) if values is None else np.asarray(values, dtype=float)
    if mask.shape != (len(qc),) or colors.shape != (len(qc),):
        raise ValueError("Spatial values and mask must have one entry per HD bin.")
    if not np.isfinite(colors).all():
        raise ValueError("Spatial HD plots require finite color values.")
    shape = (int(row.max()) + 1, int(col.max()) + 1)
    col_edges, row_edges = np.meshgrid(
        np.arange(shape[1] + 1) - 0.5, np.arange(shape[0] + 1) - 0.5
    )
    x = col_edges * transform[0, 0] + row_edges * transform[1, 0] + transform[2, 0]
    y = col_edges * transform[0, 1] + row_edges * transform[1, 1] + transform[2, 1]
    mesh_options = dict(
        shading="flat", edgecolors="none", antialiased=False, rasterized=True
    )
    background = np.full(shape, np.nan)
    if show_excluded:
        background[row[~mask], col[~mask]] = 1
    ax.pcolormesh(
        x,
        y,
        np.ma.masked_invalid(background),
        cmap=ListedColormap(["lightgray"]),
        **mesh_options,
    )
    grid = np.full(shape, np.nan)
    grid[row[mask], col[mask]] = colors[mask]
    points = ax.pcolormesh(
        x,
        y,
        np.ma.masked_invalid(grid),
        cmap=ListedColormap(["#2166ac"]) if values is None else "viridis",
        vmin=float(colors.min()),
        vmax=float(colors.max()),
        **mesh_options,
    )
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_axis_off()
    return points


def _plot_hd_qc_distributions(axes, qc: pd.DataFrame) -> None:
    """Use the same metric transforms and histograms in both QC layouts."""
    columns = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    for ax, column in zip(axes, columns):
        values = qc[column].to_numpy()
        if column in {"total_counts", "n_genes_by_counts"}:
            values = np.log10(values + 1)
            column = f"log10({column} + 1)"
        ax.hist(values, bins=80)
        ax.set_xlabel(column)
        ax.set_ylabel("Number of 8 µm bins")


def plot_hd_qc(qc: pd.DataFrame, mouse_id: str, output_dir: Path) -> None:
    """Save the notebook's distributions, spatial QC, and candidate footprints."""
    columns = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    evaluate_hd_qc(qc)
    fig, axes = plt.subplots(
        1,
        len(columns),
        figsize=(5 * len(columns), 4),
        squeeze=False,
        constrained_layout=True,
    )
    _plot_hd_qc_distributions(axes.flat, qc)
    fig.suptitle(f"{mouse_id}: QC distributions (no additional filtering)")
    _save_qc_figure(fig, output_dir / "qc_distributions.png")

    fig, axes = plt.subplots(1, 3, figsize=(21, 6), constrained_layout=True)
    for ax, column in zip(axes, columns):
        values = qc[column].to_numpy()
        label = "Mitochondrial counts (%)"
        if column != "pct_counts_mt":
            values = np.log1p(values)
            label = f"log1p({column})"
        points = plot_hd_spatial(ax, qc, np.ones(len(qc), dtype=bool), values)
        ax.set_title(label)
        fig.colorbar(points, ax=ax, shrink=0.7)
    fig.suptitle(f"{mouse_id}: spatial QC at 8 µm (no additional filtering)")
    _save_qc_figure(fig, output_dir / "spatial_qc.png")

    _, masks = evaluate_hd_qc(qc)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    for ax, (label, mask) in zip(axes.flat, masks.items()):
        points = plot_hd_spatial(ax, qc, mask, np.log1p(qc["total_counts"].to_numpy()))
        fig.colorbar(points, ax=ax, shrink=0.7, label="log1p(total_counts)")
        ax.set_title(f"{label}\n{mask.sum():,} bins ({100 * mask.mean():.1f}%)")
    fig.suptitle(f"{mouse_id}: candidate retention at 8 µm")
    _save_qc_figure(fig, output_dir / "candidate_spatial_retention.png")


def plot_hd_qc_histology(
    qc: pd.DataFrame,
    mouse_id: str,
    image: np.ndarray,
    scale: float,
    output_dir: Path,
) -> None:
    """Overlay provisional moderate QC failures on registered H&E; retain all bins."""
    evaluate_hd_qc(qc)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("The H&E scale factor must be finite and positive.")
    candidate, min_counts, min_genes, max_mt = next(
        item for item in HD_QC_CANDIDATES if item[0] == "Moderate 8 µm"
    )
    low_content = (qc["total_counts"] < min_counts) | (
        qc["n_genes_by_counts"] < min_genes
    )
    high_mt = qc["pct_counts_mt"] > max_mt
    groups = [
        ("MT-only", high_mt & ~low_content, "#d73027"),
        ("Low-content only", low_content & ~high_mt, "#0072b2"),
        ("MT + low content", high_mt & low_content, "#a000a0"),
    ]
    scaled_qc = qc.copy()
    scaled_qc[["x", "y"]] *= scale
    height, width = image.shape[:2]
    inside = scaled_qc["x"].between(-0.5, width - 0.5) & scaled_qc["y"].between(
        -0.5, height - 0.5
    )
    if not inside.any():
        raise ValueError("No bin centers overlap the H&E image; check registration.")
    fig, axes = plt.subplots(2, 2, figsize=(16, 16), constrained_layout=True)
    try:
        for ax in axes.flat:
            ax.imshow(image, origin="upper")
            ax.set_axis_off()
        axes.flat[0].set_title("H&E")
        for ax, (label, mask, color) in zip(list(axes.flat)[1:], groups):
            mesh = plot_hd_spatial(ax, scaled_qc, mask, show_excluded=False)
            mesh.set_cmap(ListedColormap([color]))
            mesh.set_alpha(0.55)
            ax.set_title(f"{label}: {mask.sum():,} bins ({100 * mask.mean():.1f}%)")
        for ax in axes.flat:
            ax.set_xlim(-0.5, width - 0.5)
            ax.set_ylim(height - 0.5, -0.5)
        fig.suptitle(
            f"{mouse_id}: {candidate} — provisional H&E review\n"
            f"Low content: UMI < {min_counts} or genes < {min_genes}; "
            f"high MT: > {max_mt:g}%\n"
            "Percentages use all input bins; no additional filtering"
        )
        _save_qc_figure(fig, output_dir / "qc_failures_he.png")
    finally:
        plt.close(fig)


def plot_hd_qc_comparison(
    candidates: pd.DataFrame,
    qc_paths: dict[str, Path],
    output_dir: Path,
) -> None:
    """Save combined distributions, retention curves, and spatial comparisons."""
    items = list(qc_paths.items())
    if items:
        fig, axes = plt.subplots(
            len(items),
            3,
            figsize=(15, 2.5 * len(items)),
            squeeze=False,
            sharex="col",
            constrained_layout=True,
        )
        for row, (sample_id, path) in zip(axes, items):
            qc = pd.read_parquet(
                path, columns=["total_counts", "n_genes_by_counts", "pct_counts_mt"]
            )
            evaluate_hd_qc(qc)
            _plot_hd_qc_distributions(row, qc)
            row[0].set_ylabel(f"{sample_id}\nNumber of 8 µm bins")
            for ax in row[1:]:
                ax.set_ylabel("")
        for ax, title in zip(axes[0], ["UMIs", "Detected genes", "Mitochondrial %"]):
            ax.set_title(title)
        for row in axes[:-1]:
            for ax in row:
                ax.set_xlabel("")
        fig.suptitle("QC distributions across mice · 8 µm · no additional filtering")
        _save_qc_figure(fig, output_dir / "cross_mouse_qc_distributions.png")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for mouse_id, frame in candidates.groupby("mouse_id", sort=False):
        for ax, column, label in zip(
            axes,
            ["bins_retained_pct", "transcripts_retained_pct"],
            ["Bins retained (%)", "Transcripts retained (%)"],
        ):
            ax.plot(
                ["None", "Lenient", "Moderate", "Higher-content"],
                frame[column],
                marker="o",
                label=mouse_id,
            )
            ax.set_ylabel(label)
    for ax in axes:
        ax.set_ylim(0, 102)
        ax.grid(alpha=0.25)
        ax.tick_params(axis="x", rotation=25)
        ax.legend(title="Mouse", fontsize="small", ncol=2)
    fig.suptitle("Candidate QC retention across 8 µm samples")
    _save_qc_figure(fig, output_dir / "cross_mouse_retention.png")

    for start in range(0, len(items), 4):
        page = items[start : start + 4]
        ncols = min(2, len(page))
        nrows = (len(page) + ncols - 1) // ncols
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(7 * ncols, 6 * nrows),
            squeeze=False,
            constrained_layout=True,
        )
        for ax, (mouse_id, path) in zip(axes.flat, page):
            qc = pd.read_parquet(path)
            _, masks = evaluate_hd_qc(qc)
            mask = masks["Moderate 8 µm"]
            plot_hd_spatial(ax, qc, mask)
            ax.set_title(f"{mouse_id}: {mask.sum():,} bins ({100 * mask.mean():.1f}%)")
        for ax in list(axes.flat)[len(page) :]:
            ax.set_visible(False)
        fig.suptitle("Moderate 8 µm: retained bins (blue), excluded bins (gray)")
        _save_qc_figure(fig, output_dir / f"cross_mouse_spatial_{start // 4 + 1}.png")


def stacked_bar(
    df: pd.DataFrame,
    value_col: str,
    title: str,
    outfile: str | Path,
    only_immune: bool = False,
) -> None:
    """Create stacked bar plot grouped by sample and condition.

    Args:
        df: DataFrame with columns 'sample', 'condition', 'celltype', and value_col.
        value_col: Column name for bar heights (e.g., 'fraction', 'count').
        title: Plot title.
        outfile: Path to save figure.
        only_immune: If True, filter to rows where celltype_is_immune=True.

    Example:
        >>> stacked_bar(df, 'fraction', 'Cell Type Distribution', 'plot.png')
    """
    d = df.copy()
    if only_immune:
        d = d[d["celltype_is_immune"]].copy()

    # Order mice by condition then sample
    order = (
        d[["sample", "condition"]]
        .drop_duplicates()
        .sort_values(
            ["condition", "sample"],
        )
    )
    d["sample"] = pd.Categorical(
        d["sample"],
        categories=order["sample"],
        ordered=True,
    )

    piv = d.pivot_table(
        index=["sample", "condition"],
        columns="celltype",
        values=value_col,
        fill_value=0.0,
        observed=True,
    ).sort_index(level=[1, 0])

    ax = piv.plot(kind="bar", stacked=True, figsize=(12, 4), width=0.9)
    ax.set_ylabel(value_col)
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    # X labels as "sample\ncondition"
    ax.set_xticklabels(
        [f"{s}\n{c}" for (s, c) in piv.index],
        rotation=0,
        ha="center",
    )

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def summary_dotplot(
    df: pd.DataFrame,
    value_col: str,
    title: str,
    outfile: str | Path,
    only_immune: bool = False,
) -> None:
    """Create dotplot with error bars (mean ± SEM) per condition and celltype.

    Args:
        df: DataFrame with columns 'condition', 'celltype', and value_col.
        value_col: Column name for plotted values.
        title: Plot title.
        outfile: Path to save figure.
        only_immune: If True, filter to immune cells only.

    Example:
        >>> summary_dotplot(df, 'fraction', 'Immune Cell Frequency', 'plot.png')
    """
    d = df.copy()
    if only_immune:
        d = d[d["celltype_is_immune"]].copy()

    # Mean ± SEM across mice, per condition and celltype
    summ = (
        d.groupby(["condition", "celltype"], observed=True)[value_col]
        .agg(mean="mean", sem=lambda x: x.std(ddof=1) / np.sqrt(len(x)))
        .reset_index()
    )

    # Simple dotplot: one panel per condition
    conditions = list(d["condition"].cat.categories)
    fig, axes = plt.subplots(
        1,
        len(conditions),
        figsize=(5 * len(conditions), 6),
        sharey=True,
    )

    if len(conditions) == 1:
        axes = [axes]

    for ax, cond in zip(axes, conditions):
        sub = summ[summ["condition"] == cond].sort_values("mean", ascending=False)
        ax.errorbar(sub["mean"], sub["celltype"], xerr=sub["sem"], fmt="o")
        ax.set_title(str(cond))
        ax.set_xlabel(value_col)
        ax.grid(True, axis="x", alpha=0.3)

    fig.suptitle(title, y=1.02)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def pretty_title(contrast: str) -> str:
    """Convert contrast key into human-readable title.

    Args:
        contrast: Contrast key formatted as '<a>_vs_<b>' (e.g., 'flash_vs_conv').

    Returns:
        Human-readable contrast title (e.g., 'FLASH vs Conventional').

    Example:
        >>> pretty_title('flash_vs_conv')
        'FLASH vs Conventional'
    """
    cond_map = {"flash": "FLASH", "conv": "Conventional", "control": "Control"}
    a, b = contrast.split("_vs_")
    return f"{cond_map.get(a, a)} vs {cond_map.get(b, b)}"
