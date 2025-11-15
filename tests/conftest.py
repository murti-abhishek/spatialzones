# tests/conftest.py
import pytest
import scanpy as sc
from importlib.resources import files
from spatialzones import datasets  # your datasets module

@pytest.fixture
def example_adata():
    """Load example spatial AnnData for tests"""
    # Use importlib.resources.files instead of pkg_resources.path
    example_file = files(datasets) / "example_spatial.h5ad"
    adata = sc.read_h5ad(example_file)
    return adata
