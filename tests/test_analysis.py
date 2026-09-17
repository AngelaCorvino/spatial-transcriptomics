"""Tests for analysis helpers."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from spatial_transcriptomics.analysis import (
    compute_statistics,
    decompose_hd_qc,
    evaluate_hd_qc,
    prepare_results,
)


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


def test_hd_candidates_use_inclusive_thresholds_and_keep_input() -> None:
    """Boundary bins match notebook retention, including transcript weighting."""
    qc = pd.DataFrame(
        {
            "total_counts": [0, 10, 25, 50, 100],
            "n_genes_by_counts": [0, 10, 20, 40, 90],
            "pct_counts_mt": [0, 25, 20, 20, 21],
        }
    )
    original = qc.copy(deep=True)
    summary, masks = evaluate_hd_qc(qc)
    summary = summary.set_index("candidate")
    assert summary["bins_retained"].tolist() == [5, 4, 2, 1]
    assert summary.loc["Moderate 8 µm", "transcripts_retained_pct"] == pytest.approx(
        100 * 75 / 185,
    )
    np.testing.assert_array_equal(
        masks["Moderate 8 µm"], [False, False, True, True, False]
    )
    pd.testing.assert_frame_equal(qc, original)


def test_hd_candidates_handle_zero_counts_and_no_retained_bins() -> None:
    qc = pd.DataFrame(
        {"total_counts": [0], "n_genes_by_counts": [0], "pct_counts_mt": [0]}
    )
    summary, _ = evaluate_hd_qc(qc)
    assert summary["transcripts_retained_pct"].tolist() == [0, 0, 0, 0]
    assert np.isnan(summary.iloc[1]["median_counts_retained"])
    with pytest.raises(ValueError, match="reference genes"):
        evaluate_hd_qc(qc.drop(columns="pct_counts_mt"))


def test_hd_failure_groups_partition_bins_and_reconcile_with_retention() -> None:
    """Exercise all eight intersections, inclusive bounds, and UMI weighting."""
    qc = pd.DataFrame(
        {
            "total_counts": [25, 24, 25, 19, 25, 24, 25, 19],
            "n_genes_by_counts": [20, 20, 19, 19, 20, 20, 19, 19],
            "pct_counts_mt": [20, 20, 20, 20, 20.1, 21, 21, 21],
        }
    )
    original = qc.copy(deep=True)
    result = decompose_hd_qc(qc)
    moderate = result[result["candidate"] == "Moderate 8 µm"]
    exclusive = moderate[moderate["group_type"] == "exclusive"]
    assert exclusive["bins"].tolist() == [1] * 8
    assert exclusive["umis"].tolist() == qc["total_counts"].tolist()
    assert exclusive["bins_pct"].sum() == pytest.approx(100)
    assert exclusive["umis_pct"].sum() == pytest.approx(100)
    marginal = moderate[moderate["group_type"] == "marginal"].set_index("group")
    assert marginal["bins"].tolist() == [4, 4, 4]
    assert marginal["umis"].tolist() == [86, 88, 93]
    candidates, _ = evaluate_hd_qc(qc)
    for _, candidate in candidates.iterrows():
        groups = result[
            (result["candidate"] == candidate["candidate"])
            & (result["group_type"] == "exclusive")
        ]
        passing = groups[groups["group"] == "passes_all"].iloc[0]
        assert groups["bins"].sum() == len(qc)
        assert groups["umis"].sum() == qc["total_counts"].sum()
        assert passing["bins"] == candidate["bins_retained"]
        assert passing["umis_pct"] == pytest.approx(
            candidate["transcripts_retained_pct"]
        )
    pd.testing.assert_frame_equal(qc, original)


def test_hd_failure_zero_umis_and_invalid_metrics() -> None:
    qc = pd.DataFrame(
        {
            "total_counts": [0],
            "n_genes_by_counts": [0],
            "pct_counts_mt": [0],
        }
    )
    result = decompose_hd_qc(qc)
    assert result["umis_pct"].isna().all()
    assert (result["umis"] == 0).all()
    moderate = result[result["candidate"] == "Moderate 8 µm"].set_index("group")
    assert moderate.loc["umi_and_genes_only", "bins"] == 1
    with pytest.raises(ValueError, match="reference genes"):
        decompose_hd_qc(qc.drop(columns="pct_counts_mt"))
    with pytest.raises(ValueError, match="finite"):
        decompose_hd_qc(qc.assign(pct_counts_mt=np.nan))
    with pytest.raises(ValueError, match="bins"):
        decompose_hd_qc(qc.iloc[:0])
