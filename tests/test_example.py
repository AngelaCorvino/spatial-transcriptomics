"""Tests for the example module."""

import pytest

from spatial_transcriptomics.example import add


def test_add():
    """Test the add function."""
    assert add(1, 2) == 3
    assert add(-1, 1) == 0