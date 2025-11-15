# spatialzones
<<<<<<< HEAD
Graph-based assignment of spatial tumor zones (inside / interface / outside) for spatial transcriptomics
=======

SpatialZones is a lightweight Python package for defining concentric spatial regions
(**inside → interface → outside**) around tumor cells using graph distance on
spatial transcriptomics data. It integrates seamlessly with **Scanpy**, **AnnData**, and **Squidpy**.

The package provides:
- `assign_tumor_regions()` — core algorithm for region assignment  
- `plot_region_gene_spatial()` — spatial gene + region visualization  
- `plot_region_composition()` — composition barplots  
- `plot_gene_violin_by_region()` — violin plots  
- `plot_graph_vs_euclidean()` — correlation between Euclidean and graph distance  

This package is designed for spatial transcriptomics applications
(Visium, Xenium, CosMx, and custom ST assays).

---

## 📦 Installation

```bash
pip install spatialzones
>>>>>>> ac84a1f1 (Initial commit)
