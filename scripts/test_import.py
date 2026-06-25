#!/usr/bin/env python3
"""Check that the cluster environment can import project dependencies."""

from __future__ import annotations

import sys

import anndata
import matplotlib
import pandas
import scanpy
import squidpy
import yaml

import spatial_transcriptomics


def main() -> int:
    """Print import and version information for a minimal cluster smoke test."""
    python_version = sys.version.replace("\n", " ")

    print(f"Python version: {python_version}")
    print(f"scanpy version: {scanpy.__version__}")
    print(f"anndata version: {anndata.__version__}")
    print(f"pandas version: {pandas.__version__}")
    print(f"matplotlib version: {matplotlib.__version__}")
    print(f"PyYAML version: {yaml.__version__}")
    print(f"squidpy version: {squidpy.__version__}")
    print(f"spatial_transcriptomics version: {spatial_transcriptomics.__version__}")
    print("spatial_transcriptomics import: OK")
    print("Cluster test completed successfully")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
