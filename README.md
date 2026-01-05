# spatialzones

[![Python Package](https://github.com/murti-abhishek/spatialzones/actions/workflows/python-package.yml/badge.svg)](https://github.com/murti-abhishek/spatialzones/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/murti-abhishek/spatialzones/branch/main/graph/badge.svg?token=YOUR_CODECOV_TOKEN)](https://codecov.io/gh/murti-abhishek/spatialzones)

**`spatialzones`** is a lightweight Python package for assigning and analyzing **tumor-centric spatial zones** in spatial transcriptomics data.

---

## Overview

Understanding the tumor microenvironment (TME) requires precise spatial context. `spatialzones` leverages graph-based distances through cell-cell connectivity networks to identify:

- *Inside tumor*: Core tumor regions and immediately adjacent cells
- *Interface*: Cells at the tumor boundary with direct tumor contact
- *Outside tumor*: Cells in the surrounding stroma without tumor proximity

Unlike simple Euclidean distance approaches, graph-based methods respect tissue connectivity and provide biologically meaningful region assignments.

`spatialzones` provides a **deterministic, annotation-aware framework** to:

- define tumor-centric spatial regions,
- retain interpretability across samples and cohorts

The package was developed with **Xenium** and **Visium-like** spatial data in mind, but operates on standard `AnnData` objects.

---

## Key Features

- Graph-based region assignment: Uses Dijkstra's algorithm on spatial connectivity graphs
- Dual distance metrics: Computes both graph distance and Euclidean distance for comparison
- Flexible hop parameters: Customizable thresholds for defining region boundaries
- Rich visualizations: Built-in plotting for spatial distributions, compositions, and gene expression
- Scanpy/Squidpy integration: Works seamlessly with AnnData objects
- Statistical comparisons: Compare observed vs expected regional compositions

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/murti-abhishek/spatialzones.git
cd spatialzones
pip install -e .
```

### Dependencies

Core dependencies include:

- `scanpy`
- `squidpy`
- `numpy`, `pandas`
- `scipy`
- `matplotlib`, `seaborn`

---

## Quick start

```python
import scanpy as sc
from spatialzones.assign import assign_tumor_regions
from spatialzones.plot import plot_region_gene_spatial
from spatialzones.plot import plot_expected_vs_observed_regions
```

### Load example dataset

```python
adata = sc.read_h5ad("src/spatialzones/datasets/example_spatial.h5ad")
```

### Assign tumor regions

```python
assign_tumor_regions(
    adata,
    tumor_type="Tumor Hepatoblast", ## any cell type in your annotation_key
    annotation_key="temp_annotations" ## metadata that stores your annotations
)
```

### Plot regions and gene expression

```python
plot_region_gene_spatial(adata, gene="AHSG") ## AHSG is expressed by Hepatoblast like tumor cells
```

### Plot expected vs observed celltype distribution

```python
sz.plot_expected_vs_observed_regions(
    adata,
    remove_celltypes=["Tumor Hepatoblast"] ## you can visualize the celltype distribution with or without the tumor cells
)
```
---

## What `assign_tumor_regions` adds

After running `assign_tumor_regions`, the following fields are populated:

### `adata.obs`

| Column | Description |
|------|-------------|
| `region` | Categorical region assignment (`inside`, `interface`, `outside`) |
| `graph_dist_to_tumor` | Graph distance to nearest tumor cell |
| `euclid_dist_to_tumor` | Euclidean distance to nearest tumor cell |

### `adata.uns`

| Key | Description |
|----|-------------|
| `region_colors` | Consistent color mapping for regions |

These outputs are intentionally simple and interoperable with Scanpy/Squidpy workflows.

---

## Plotting utilities

Plotting functions live in `spatialzones.plot` and include:

- Spatial visualization of regions and genes
- Gene expression by region (violin, dotplot)
- Graph vs Euclidean distance comparisons
- Plotting expected vs actual celltype distribution to identify cell types enriched in specific regions

Plotting is **optional** and kept separate from core logic.

---

## When to use spatialzones

`spatialzones` is well-suited for:

- Tumor microenvironment analyses
- Comparing spatial organization across samples
- Integrating spatial structure into downstream modeling

---

## Citation

If you use `spatialzones` in your work, please cite or acknowledge the repository.  
A methods-style description can be provided upon request.

---

## Author

**Abhishek Murti**  
GitHub: [@murti-abhishek](https://github.com/murti-abhishek)