"""Configuration management for spatial transcriptomics analysis.

Loads configuration from local_config.yaml or provides sensible defaults.
Avoids hardcoded absolute paths in notebooks.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:
    yaml = None


def load_config() -> dict[str, Any]:
    """Load configuration from local_config.yaml or return defaults.

    Returns:
        Configuration dictionary with data paths and analysis settings.
    """
    repo_root = Path(__file__).parent.parent.parent
    local_config_path = repo_root / "local_config.yaml"

    config: dict[str, Any] = {
        "repo_root": str(repo_root),
        "data_dir": str(repo_root / "data"),
        "output_dir": str(repo_root / "results"),
        "cache_dir": str(repo_root / ".cache"),
    }

    if local_config_path.exists() and yaml is not None:
        try:
            with local_config_path.open() as f:
                local_config = yaml.safe_load(f) or {}
                config.update(local_config)
        except (OSError, yaml.YAMLError) as exc:
            warnings.warn(
                f"Failed to load local_config.yaml: {exc}",
                stacklevel=2,
            )

    return config


def get_data_path(filename: str) -> Path:
    """Get absolute path for a data file.

    Args:
        filename: Name of the data file.

    Returns:
        Absolute path to the data file.
    """
    config = load_config()
    return Path(config["data_dir"]) / filename


def get_output_path(filename: str) -> Path:
    """Get absolute path for an output file.

    Args:
        filename: Name of the output file.

    Returns:
        Absolute path to the output file.
    """
    config = load_config()
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / filename
