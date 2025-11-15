# tests/test_plot.py
import pytest
import matplotlib.pyplot as plt
from spatialzones import plot
from spatialzones.assign import assign_tumor_regions  # assuming you have this

def test_plot_region_gene_spatial_basic(example_adata):
    adata = example_adata

    # Pick a tumor type that exists
    tumor_type = adata.obs["temp_annotations"].unique()[0]

    # Assign regions for plotting
    adata = assign_tumor_regions(
        adata,
        tma_cores=[adata.obs["tma_core"].unique()[0]],
        tumor_type=tumor_type,
    )

    # Test spatial gene plotting
    gene = adata.var_names[0]  # pick first gene
    plot.plot_region_gene_spatial(adata, gene)  # should not raise

def test_plot_region_composition_basic(example_adata):
    adata = example_adata

    # Assign regions for plotting
    tumor_type = adata.obs["temp_annotations"].unique()[0]
    adata = assign_tumor_regions(
        adata,
        tma_cores=[adata.obs["tma_core"].unique()[0]],
        tumor_type=tumor_type,
    )

    # Test composition plotting
    plot.plot_region_composition(adata)  # should not raise

def test_plot_gene_violin_by_region_basic(example_adata):
    adata = example_adata
    tumor_type = adata.obs["temp_annotations"].unique()[0]

    adata = assign_tumor_regions(
        adata,
        tma_cores=[adata.obs["tma_core"].unique()[0]],
        tumor_type=tumor_type,
    )

    gene = adata.var_names[0]
    plot.plot_gene_violin_by_region(adata, gene)

def test_plot_gene_dotplot_by_region_basic(example_adata):
    adata = example_adata
    tumor_type = adata.obs["temp_annotations"].unique()[0]

    adata = assign_tumor_regions(
        adata,
        tma_cores=[adata.obs["tma_core"].unique()[0]],
        tumor_type=tumor_type,
    )

    gene = adata.var_names[0]
    plot.plot_gene_dotplot_by_region(adata, gene)
