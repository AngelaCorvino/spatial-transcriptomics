"""Tests for analysis helpers."""

import argparse
import hashlib
import json
import runpy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from spatial_transcriptomics.analysis import (
    compute_statistics,
    decompose_hd_qc,
    evaluate_hd_qc,
    prepare_results,
    select_hd_qc_bins,
    selected_hd_qc_policy,
    summarize_annotation_markers,
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


def test_selected_hd_filter_preserves_boundaries_and_has_no_upper_cutoffs() -> None:
    qc = pd.DataFrame(
        {
            "total_counts": [25, 24, 25, 25, 100000],
            "n_genes_by_counts": [20, 20, 19, 20, 7000],
            "pct_counts_mt": [20, 20, 20, 20.1, 0],
            "pct_counts_rp": [99] * 5,
        },
        index=["boundary", "umi", "genes", "mt", "high_content"],
    )
    original = qc.copy(deep=True)
    selected = select_hd_qc_bins(qc)
    assert selected["qc_keep"].tolist() == [True, False, False, False, True]
    assert selected["fails_umi"].tolist() == [False, True, False, False, False]
    assert selected["fails_genes"].tolist() == [False, False, True, False, False]
    assert selected["fails_mt"].tolist() == [False, False, False, True, False]
    pd.testing.assert_frame_equal(qc, original)
    assert selected_hd_qc_policy()["anatomical_roi_applied"] is False
    with pytest.raises(ValueError, match="unique"):
        select_hd_qc_bins(pd.concat([qc, qc.iloc[:1]]))
    with pytest.raises(ValueError, match="negative"):
        select_hd_qc_bins(qc.assign(total_counts=-1))


def _annotation_fixture(mouse):
    """Synthetic three-region tissue with raw counts and an affine spatial grid."""
    ad = pytest.importorskip("anndata")
    from scipy import sparse

    from spatial_transcriptomics.data import apply_hd_qc_filter

    rng = np.random.default_rng(42)
    counts = rng.poisson(1.5, size=(90, 64)).astype(np.float32)
    counts[:, 0] = 1
    for group in range(3):
        counts[group * 30 : (group + 1) * 30, 1 + group * 10 : 11 + group * 10] += 8
    counts[:, -3:] = 0
    counts[:30, -3] = 12  # Cd3e
    counts[30:60, -2] = 12  # Acta1
    row, col = np.divmod(np.arange(90), 10)
    names = [f"s_008um_{r:05d}_{c:05d}-1" for r, c in zip(row, col)]
    genes = ["mt-Nd1"] + [f"Gene{i}" for i in range(60)] + ["Cd3e", "Acta1", "Zero"]
    adata = ad.AnnData(
        sparse.csr_matrix(counts),
        obs=pd.DataFrame({"in_tissue": 1, "sample_id": mouse}, index=names),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm["spatial"] = np.column_stack([100 + 8 * col, 100 + 8 * row]).astype(
        float
    )
    adata.uns["visium_hd"] = {"bin_size_um": 8, "source_path": "synthetic.tar.gz"}
    return apply_hd_qc_filter(adata, ["mt-Nd1"])


@pytest.mark.parametrize("resolutions", [[0.2], [0.1, 0.2, 0.3]])
def test_hd_annotation_run_preserves_counts_and_writes_review_evidence(
    tmp_path, monkeypatch, resolutions
):
    """Exercise real sparse normalization, clustering, marker ranking and plots."""
    pytest.importorskip("igraph")
    import anndata as ad

    repo = Path(__file__).resolve().parents[1]
    script = runpy.run_path(str(repo / "scripts/03_integration.py"))
    input_dir = tmp_path / "filtered"
    output_dir = tmp_path / "annotation"
    he_root = tmp_path / "histology"
    originals = {}
    summary = []
    hashes = {}
    for mouse in ["FD1", "FD2"]:
        adata, selection = _annotation_fixture(mouse)
        originals[mouse] = adata.copy()
        sample_dir = input_dir / mouse
        sample_dir.mkdir(parents=True)
        adata.write_h5ad(sample_dir / "filtered_counts.h5ad")
        selection.iloc[::-1].to_parquet(sample_dir / "bin_selection.parquet")
        hashes[mouse] = hashlib.sha256(
            (sample_dir / "filtered_counts.h5ad").read_bytes()
        ).hexdigest()
        (he_root / mouse).mkdir(parents=True)
        (he_root / mouse / "spatial.tar.gz").touch()
        summary.append(
            {
                "mouse_id": mouse,
                "bins_retained": adata.n_obs,
                "features_retained": adata.n_vars,
            }
        )
    (input_dir / "run.json").write_text(
        json.dumps({"status": "complete", "completed_mice": ["FD1", "FD2"]})
    )
    pd.DataFrame(summary).to_csv(input_dir / "filter_summary.csv", index=False)
    marker_path = tmp_path / "markers.csv"
    marker_table = pd.DataFrame(
        [
            ("T", "Cd3e", "synthetic"),
            ("Muscle", "Acta1", "synthetic"),
            ("Coverage", "Zero", "synthetic"),
            ("Coverage", "Absent", "synthetic"),
            ("Shared", "Cd3e", "synthetic"),
        ],
        columns=["panel", "gene", "source"],
    )
    marker_table.to_csv(marker_path, index=False)
    args = argparse.Namespace(
        all_mice=False,
        mice=["FD1", "FD2"],
        input_dir=input_dir,
        output_dir=output_dir,
        markers=marker_path,
        he_source_root=he_root,
        n_top_genes=40,
        n_pcs=5,
        n_neighbors=10,
        resolution=None,
        resolutions=resolutions,
        seed=0,
    )
    config = {"repo_root": str(repo), "data_dir": str(tmp_path / "raw")}
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    run = script["run_visium_hd_annotation"]
    with pytest.raises(ValueError, match="SLURM"):
        run(args, config)
    monkeypatch.setenv("SLURM_JOB_ID", "synthetic-test")
    monkeypatch.setitem(
        run.__globals__,
        "load_visium_hd_histology",
        lambda *_: (np.ones((220, 220, 3)), 1.0),
    )
    embeddings = {}
    cluster = script["cluster_and_umap"]

    def capture_embedding(adata, **kwargs):
        mouse = adata.obs["sample_id"].iloc[0]
        assert mouse not in embeddings  # One PCA/graph/UMAP per mouse.
        cluster(adata, **kwargs)
        embeddings[mouse] = (
            adata.obsm["X_pca"].copy(),
            adata.obsm["X_umap"].copy(),
            adata.obsp["connectivities"].copy(),
        )

    monkeypatch.setitem(run.__globals__, "cluster_and_umap", capture_embedding)
    assert run(args, config) == 0
    manifest = json.loads((output_dir / "run.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["completed_mice"] == ["FD1", "FD2"]
    assert manifest["batch_correction"] is False
    assert manifest["settings"]["resolutions"] == resolutions
    assert manifest["shared_embedding"] is True
    assert manifest["markers"]["table"] == marker_table.to_dict("records")
    assert (
        manifest["markers"]["sha256"]
        == hashlib.sha256(marker_path.read_bytes()).hexdigest()
    )
    summary = pd.read_csv(output_dir / "annotation_summary.csv")
    assert summary["marker_genes_present"].eq(3).all()
    assert len(summary) == 2 * len(resolutions)
    assert set(summary["resolution"]) == set(resolutions)
    for mouse in ["FD1", "FD2"]:
        sample = output_dir / mouse
        restored = ad.read_h5ad(sample / "clustered_bins.h5ad")
        original = originals[mouse]
        assert (restored.layers["counts"] != original.X).nnz == 0
        assert restored.obs_names.equals(original.obs_names)
        assert restored.var_names.equals(original.var_names)
        np.testing.assert_array_equal(
            restored.obsm["spatial"], original.obsm["spatial"]
        )
        expected = original.X.toarray()
        expected = np.log2(1 + expected * 10000 / expected.sum(axis=1)[:, None])
        np.testing.assert_allclose(restored.X.toarray(), expected, rtol=1e-6)
        assert not restored.var.loc["mt-Nd1", "highly_variable"]
        pca, umap, graph = embeddings[mouse]
        np.testing.assert_array_equal(restored.obsm["X_pca"], pca)
        np.testing.assert_array_equal(restored.obsm["X_umap"], umap)
        assert (restored.obsp["connectivities"] != graph).nnz == 0
        assignments = pd.read_csv(
            sample / "cluster_assignments.csv", index_col="barcode", dtype=str
        )
        assert assignments.index.equals(original.obs_names)
        assert assignments.columns.tolist() == [f"leiden_r{r:g}" for r in resolutions]
        if len(resolutions) > 1:
            assert "leiden" not in restored.obs
            assert (sample / "resolution_comparison.png").stat().st_size > 1000
        else:
            assert restored.obs["leiden"].nunique() >= 2
        assert (
            restored.uns["tissue_annotation"]["status"] == "clusters_ready_for_review"
        )
        assert (
            hashlib.sha256(
                (input_dir / mouse / "filtered_counts.h5ad").read_bytes()
            ).hexdigest()
            == hashes[mouse]
        )
        coverage = pd.read_csv(sample / "marker_coverage.csv").set_index("gene")
        assert not coverage.loc["Absent", "in_matrix"]
        assert pd.isna(coverage.loc["Absent", "bins_expressing"])
        assert coverage.loc["Zero", "in_matrix"]
        assert coverage.loc["Zero", "bins_expressing"] == 0
        for resolution in resolutions:
            key = f"leiden_r{resolution:g}"
            assert restored.obs[key].nunique() >= 2
            assert assignments[key].tolist() == restored.obs[key].astype(str).tolist()
            review = sample
            if len(resolutions) > 1:
                review = sample / f"resolution_{resolution:g}"
            labels = pd.read_csv(review / "cluster_annotation.csv")
            assert labels["bins"].sum() == original.n_obs
            assert labels["resolution"].eq(resolution).all()
            assert labels["annotation_status"].eq("unreviewed").all()
            assert labels["tissue_annotation"].isna().all()
            sizes = restored.obs[key].value_counts()
            row = summary.loc[
                summary["mouse_id"].eq(mouse) & summary["resolution"].eq(resolution)
            ].iloc[0]
            assert row["clusters"] == len(sizes)
            assert row["bins_in_clusters_below_100"] == sizes[sizes < 100].sum()
            for name in ["cluster_markers.csv", "marker_expression.csv"]:
                assert (review / name).stat().st_size > 0
            for name in ["clusters.png", "marker_dotplot.png"]:
                assert (review / name).stat().st_size > 1000
        for name in ["marker_spatial.png", "qc_failures_he.png"]:
            assert (sample / name).stat().st_size > 1000
    with pytest.raises(FileExistsError):
        run(args, config)
    args.output_dir = tmp_path / "failed"
    args.mice = ["FD1"]
    changed = originals["FD1"].copy()
    changed.X.data /= 2
    changed.write_h5ad(input_dir / "FD1" / "filtered_counts.h5ad")
    with pytest.raises(ValueError, match="integer counts"):
        run(args, config)
    failed = json.loads((args.output_dir / "run.json").read_text())
    assert failed["status"] == "failed" and not failed["completed_mice"]
    assert not (args.output_dir / "FD1" / "clustered_bins.h5ad").exists()


def test_hd_annotation_defaults_to_requested_resolution_comparison(monkeypatch):
    repo = Path(__file__).resolve().parents[1]
    script = runpy.run_path(str(repo / "scripts/03_integration.py"))
    monkeypatch.setattr("sys.argv", ["03_integration.py", "--visium-hd"])
    args = script["parse_args"]()
    assert args.resolution is None
    assert args.resolutions == [0.1, 0.2, 0.3]


@pytest.mark.parametrize(
    "contents, message",
    [
        ("panel,gene\nT,Cd3e\n", "columns"),
        ("panel,gene,source\n", "no empty fields"),
        ("panel,gene,source\nT, ,reference\n", "no empty fields"),
        ("panel,gene,source\nT,Cd3e,\n", "no empty fields"),
        (
            "panel,gene,source\nT,Cd3e,reference\nT,Cd3e,other reference\n",
            "unique within",
        ),
    ],
)
def test_annotation_marker_csv_rejects_incomplete_or_duplicate_entries(
    tmp_path, contents, message
):
    repo = Path(__file__).resolve().parents[1]
    script = runpy.run_path(str(repo / "scripts/03_integration.py"))
    marker_path = tmp_path / "markers.csv"
    marker_path.write_text(contents)
    with pytest.raises(ValueError, match=message):
        script["load_annotation_markers"](marker_path)


def test_annotation_marker_expression_uses_unscaled_values_and_detection():
    import anndata as ad
    from scipy import sparse

    data = ad.AnnData(
        sparse.csr_matrix([[0.0, 2.0], [2.0, 4.0], [6.0, 0.0], [0.0, 0.0]]),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["0", "0", "1", "1"])}, index=["a", "b", "c", "d"]
        ),
        var=pd.DataFrame(index=["A", "B"]),
    )
    coverage, expression = summarize_annotation_markers(
        data, {"panel": ["A", "B", "missing"]}
    )
    expression = expression.set_index(["cluster", "gene"])
    assert expression.loc[("0", "A"), "mean_log_expression"] == 1
    assert expression.loc[("1", "A"), "mean_log_expression"] == 3
    assert expression.loc[("0", "A"), "bins_expressing_pct"] == 50
    assert expression.loc[("1", "B"), "bins_expressing_pct"] == 0
    assert not coverage.set_index("gene").loc["missing", "in_matrix"]
