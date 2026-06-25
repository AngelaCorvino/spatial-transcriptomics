"""Tests for configuration helpers."""

from pathlib import Path

import pytest

from spatial_transcriptomics import config as config_module


def test_load_config_returns_repo_defaults(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Default config should be derived from the repository layout."""
    repo_root = tmp_path / "repo"
    module_path = repo_root / "src" / "spatial_transcriptomics" / "config.py"
    module_path.parent.mkdir(parents=True)
    module_path.touch()
    monkeypatch.setattr(config_module, "__file__", str(module_path))

    config = config_module.load_config()

    assert config == {
        "repo_root": str(repo_root),
        "data_dir": str(repo_root / "data"),
        "output_dir": str(repo_root / "results"),
        "cache_dir": str(repo_root / ".cache"),
    }


def test_load_config_merges_local_yaml(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Local YAML should override defaults and allow nested settings."""
    repo_root = tmp_path / "repo"
    custom_data_dir = tmp_path / "custom-data"
    module_path = repo_root / "src" / "spatial_transcriptomics" / "config.py"
    module_path.parent.mkdir(parents=True)
    module_path.touch()
    (repo_root / "local_config.yaml").write_text(
        f"data_dir: {custom_data_dir}\nanalysis:\n  resolution: 1.2\n",
    )
    monkeypatch.setattr(config_module, "__file__", str(module_path))

    config = config_module.load_config()

    assert config["data_dir"] == str(custom_data_dir)
    assert config["analysis"] == {"resolution": 1.2}


def test_load_config_accepts_explicit_config_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Explicit config files should be resolved from the repository root."""
    repo_root = tmp_path / "repo"
    module_path = repo_root / "src" / "spatial_transcriptomics" / "config.py"
    module_path.parent.mkdir(parents=True)
    module_path.touch()
    config_dir = repo_root / "configs"
    config_dir.mkdir()
    (config_dir / "local.yaml").write_text("output_dir: outputs\n")
    monkeypatch.setattr(config_module, "__file__", str(module_path))

    config = config_module.load_config("configs/local.yaml")

    assert config["output_dir"] == str(repo_root / "outputs")


def test_resolve_config_path_uses_repo_root(tmp_path: Path) -> None:
    """Relative config paths should resolve from repo_root."""
    config = {"repo_root": str(tmp_path)}

    result = config_module.resolve_config_path(config, "results/qc")

    assert result == tmp_path / "results" / "qc"


def test_get_config_section_rejects_non_dictionary() -> None:
    """Nested sections should fail clearly when YAML has the wrong shape."""
    with pytest.raises(TypeError, match="paths"):
        config_module.get_config_section({"paths": "results"}, "paths")


def test_get_output_path_creates_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Output helper should create the target directory on demand."""
    output_dir = tmp_path / "results"
    monkeypatch.setattr(
        config_module,
        "load_config",
        lambda: {
            "repo_root": str(tmp_path),
            "data_dir": str(tmp_path / "data"),
            "output_dir": str(output_dir),
            "cache_dir": str(tmp_path / ".cache"),
        },
    )

    result = config_module.get_output_path("plot.png")

    assert result == output_dir / "plot.png"
    assert output_dir.exists()


def test_get_data_path_rejects_non_string_config(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Path helpers should reject invalid config types early."""
    monkeypatch.setattr(
        config_module,
        "load_config",
        lambda: {
            "repo_root": str(tmp_path),
            "data_dir": {"wrong": "type"},
            "output_dir": str(tmp_path / "results"),
            "cache_dir": str(tmp_path / ".cache"),
        },
    )

    with pytest.raises(TypeError, match="data_dir"):
        config_module.get_data_path("input.h5ad")
