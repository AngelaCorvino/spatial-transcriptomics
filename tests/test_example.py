"""Tests for the example module."""


from spatial_transcriptomics.example import add


def test_add() -> None:
    """Test the add function."""
    assert add(1, 2) == 3
    assert add(-1, 1) == 0
