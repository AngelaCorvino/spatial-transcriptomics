"""Tests for plotting helpers."""

from pathlib import Path

import matplotlib as mpl
import pandas as pd

from spatial_transcriptomics.plotting import pretty_title, stacked_bar, summary_dotplot

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
