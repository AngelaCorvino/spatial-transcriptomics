"""Plotting utilities for spatial transcriptomics analysis.

Provides functions for visualization of cell type composition,
differential expression results, and other analysis outputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from spatial_transcriptomics.analysis import evaluate_hd_qc


def _save_qc_figure(fig, outfile: Path) -> None:
    """Save and immediately release a batch QC figure."""
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _qc_spatial_axis(ax, qc: pd.DataFrame, mask, values=None):
    """Draw all bins, with excluded bins gray and image-style orientation."""
    ax.scatter(
        qc.loc[~mask, "x"],
        qc.loc[~mask, "y"],
        c="lightgray",
        s=0.3,
        linewidths=0,
        rasterized=True,
    )
    points = ax.scatter(
        qc.loc[mask, "x"],
        qc.loc[mask, "y"],
        c="#2166ac" if values is None else values[mask],
        **({} if values is None else {"cmap": "viridis"}),
        s=0.3,
        linewidths=0,
        rasterized=True,
    )
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_axis_off()
    return points


def plot_hd_qc(qc: pd.DataFrame, mouse_id: str, output_dir: Path) -> None:
    """Save the notebook's distributions, spatial QC, and candidate footprints."""
    columns = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_rp"]
    columns = [column for column in columns if column in qc]
    fig, axes = plt.subplots(
        1,
        len(columns),
        figsize=(5 * len(columns), 4),
        squeeze=False,
        constrained_layout=True,
    )
    for ax, column in zip(axes.flat, columns):
        values = qc[column].to_numpy()
        if column in {"total_counts", "n_genes_by_counts"}:
            values = np.log10(values + 1)
            column = f"log10({column} + 1)"
        ax.hist(values, bins=80)
        ax.set_xlabel(column)
        ax.set_ylabel("Number of 8 µm bins")
    fig.suptitle(f"{mouse_id}: QC distributions (no additional filtering)")
    _save_qc_figure(fig, output_dir / "qc_distributions.png")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    for ax, column in zip(axes, columns[:2]):
        points = _qc_spatial_axis(
            ax, qc, np.ones(len(qc), dtype=bool), np.log1p(qc[column].to_numpy())
        )
        ax.set_title(f"log1p({column})")
        fig.colorbar(points, ax=ax, shrink=0.7)
    fig.suptitle(f"{mouse_id}: spatial QC at 8 µm")
    _save_qc_figure(fig, output_dir / "spatial_qc.png")

    _, masks = evaluate_hd_qc(qc)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    for ax, (label, mask) in zip(axes.flat, masks.items()):
        points = _qc_spatial_axis(ax, qc, mask, np.log1p(qc["total_counts"].to_numpy()))
        fig.colorbar(points, ax=ax, shrink=0.7, label="log1p(total_counts)")
        ax.set_title(f"{label}\n{mask.sum():,} bins ({100 * mask.mean():.1f}%)")
    fig.suptitle(f"{mouse_id}: candidate retention at 8 µm")
    _save_qc_figure(fig, output_dir / "candidate_spatial_retention.png")


def plot_hd_qc_comparison(
    candidates: pd.DataFrame,
    qc_paths: dict[str, Path],
    output_dir: Path,
) -> None:
    """Save cross-mouse retention curves and a paginated moderate-filter grid."""
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

    items = list(qc_paths.items())
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
            _qc_spatial_axis(ax, qc, mask)
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
