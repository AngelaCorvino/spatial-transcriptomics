"""Tests for analysis helpers."""

from types import SimpleNamespace

import pandas as pd

from spatial_transcriptomics.analysis import compute_statistics, prepare_results


def test_compute_statistics_returns_expected_summary() -> None:
    """Statistics should reflect observation, feature, and cluster counts."""
    adata = SimpleNamespace(
        n_obs=3,
        n_vars=4,
        obs=pd.DataFrame(
            {
                "total_counts": [10, 20, 30],
                "leiden": ["0", "1", "1"],
            },
        ),
    )

    stats = compute_statistics(adata)

    assert stats == {
        "n_obs": 3,
        "n_vars": 4,
        "mean_counts": 20.0,
        "n_clusters": 2,
    }


def test_prepare_results_embeds_summary_and_payload() -> None:
    """Prepared results should preserve the payload and add a human summary."""
    adata = SimpleNamespace(name="toy")
    stats = {"n_obs": 5, "n_vars": 2, "mean_counts": 4.0}

    result = prepare_results(adata, stats)

    assert result["data"] is adata
    assert result["statistics"] == stats
    assert result["summary"] == "Analysis complete for 5 observations and 2 genes."
