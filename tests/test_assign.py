# tests/test_assign.py
import pytest
from spatialzones.assign import assign_tumor_regions  # your function
from spatialzones.plot import plot_region_composition  # optional: for additional checks

def test_assign_basic(example_adata):
    adata = example_adata

    # Pick a tumor type that exists
    tumor_type = adata.obs["temp_annotations"].unique()[0]

    # Assign regions
    adata = assign_tumor_regions(
        adata,
        tma_cores=[adata.obs["tma_core"].unique()[0]],
        tumor_type=tumor_type,
    )

    # Check that regions were assigned
    assert 'region' in adata.obs, "Region column not found after assignment."
    assert 'region_colors' in adata.uns, "Region colors not set in .uns after assignment."

    # Optional: quick check that plot works without error
    # This avoids needing a separate test if you want light integration testing
    plot_region_composition(adata)
