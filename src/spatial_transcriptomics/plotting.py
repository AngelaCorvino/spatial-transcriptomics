"""Plotting utilities for spatial transcriptomics analysis.

Provides functions for visualization of cell type composition,
differential expression results, and other analysis outputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    order = d[["sample", "condition"]].drop_duplicates().sort_values(
        ["condition", "sample"],
    )
    d["sample"] = pd.Categorical(
        d["sample"],
        categories=order["sample"],
        ordered=True,
    )

    piv = (
        d.pivot_table(
            index=["sample", "condition"],
            columns="celltype",
            values=value_col,
            fill_value=0.0,
            observed=True,
        )
        .sort_index(level=[1, 0])
    )

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

    for ax, cond in zip(axes, conditions, strict=True):
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
