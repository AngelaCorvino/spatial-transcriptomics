"""Configuration management for spatial transcriptomics analysis.

Loads configuration from YAML files and provides sensible defaults.
Avoids hardcoded absolute paths in notebooks and scripts.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from pathlib import Path

try:
    import yaml as yaml_module
except ImportError:
    yaml_module = None  # type: ignore[assignment]


def _resolve_path(repo_root: Path, path: str) -> str:
    path_obj = Path(path)
    if path_obj.is_absolute():
        return str(path_obj)
    return str((repo_root / path_obj).resolve())


def resolve_config_path(config: Mapping[str, object], raw_path: str | Path) -> Path:
    """Resolve a path from config relative to the repository root."""
    path = Path(raw_path)
    if path.is_absolute():
        return path

    repo_root = config.get("repo_root", ".")
    if not isinstance(repo_root, (str, Path)):
        raise TypeError("Configuration value 'repo_root' must be a string path.")

    return (Path(repo_root) / path).resolve()


def get_config_section(
    config: Mapping[str, object],
    key: str,
) -> dict[str, object]:
    """Return a nested config section and validate it is a dictionary."""
    value = config.get(key, {})
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise TypeError(f"Configuration value '{key}' must be a dictionary.")
    return {
        str(section_key): section_value
        for section_key, section_value in value.items()
    }


def get_config_path(
    config: Mapping[str, object],
    key: str,
    default: str | Path,
) -> Path:
    """Return a top-level path config entry resolved to an absolute path."""
    value = config.get(key, default)
    if not isinstance(value, (str, Path)):
        raise TypeError(f"Configuration value '{key}' must be a string path.")
    return resolve_config_path(config, value)


def get_section_path(
    config: Mapping[str, object],
    section: Mapping[str, object],
    key: str,
    default: str | Path,
) -> Path:
    """Return a nested path config entry resolved to an absolute path."""
    value = section.get(key, default)
    if not isinstance(value, (str, Path)):
        raise TypeError(f"Configuration value '{key}' must be a string path.")
    return resolve_config_path(config, value)


def load_config(config_path: str | Path | None = None) -> dict[str, object]:
    """Load configuration from YAML file or return defaults.

    Args:
        config_path: Optional path to a config YAML file. If omitted, uses
            `local_config.yaml` in the repository root.

    Returns:
        Configuration dictionary with data paths and analysis settings.
    """
    repo_root = Path(__file__).resolve().parent.parent.parent
    default_config: dict[str, object] = {
        "repo_root": str(repo_root),
        "data_dir": str(repo_root / "data"),
        "output_dir": str(repo_root / "results"),
        "cache_dir": str(repo_root / ".cache"),
    }

    if config_path is not None:
        config_file = Path(config_path)
        if not config_file.is_absolute():
            config_file = repo_root / config_file
    else:
        config_file = repo_root / "local_config.yaml"

    if config_file.exists():
        if yaml_module is None:
            raise RuntimeError("PyYAML is required to read configuration files.")
        try:
            with config_file.open() as f:
                loaded_config = yaml_module.safe_load(f) or {}
        except (OSError, yaml_module.YAMLError) as exc:
            warnings.warn(
                f"Failed to load config file {config_file}: {exc}",
                stacklevel=2,
            )
            loaded_config = {}
        default_config.update(loaded_config)
    elif config_path is not None:
        raise FileNotFoundError(f"Config file not found: {config_file}")

    for key in ["data_dir", "output_dir", "cache_dir"]:
        value = default_config.get(key)
        if isinstance(value, str):
            default_config[key] = _resolve_path(repo_root, value)

    return default_config


def _require_path_value(config: dict[str, object], key: str) -> str:
    """Return a path config entry and validate it is a string."""
    value = config[key]
    if not isinstance(value, str):
        message = f"Configuration value '{key}' must be a string path."
        raise TypeError(message)
    return value


def get_data_path(filename: str) -> Path:
    """Get absolute path for a data file.

    Args:
        filename: Name of the data file.

    Returns:
        Absolute path to the data file.
    """
    config = load_config()
    return Path(_require_path_value(config, "data_dir")) / filename


def get_output_path(filename: str) -> Path:
    """Get absolute path for an output file.

    Args:
        filename: Name of the output file.

    Returns:
        Absolute path to the output file.
    """
    config = load_config()
    output_dir = Path(_require_path_value(config, "output_dir"))
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / filename
