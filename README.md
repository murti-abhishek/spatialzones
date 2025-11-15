# spatialzones

[![Python Package](https://github.com/murti-abhishek/spatialzones/actions/workflows/python-package.yml/badge.svg)](https://github.com/murti-abhishek/spatialzones/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/murti-abhishek/spatialzones/branch/main/graph/badge.svg?token=YOUR_CODECOV_TOKEN)](https://codecov.io/gh/murti-abhishek/spatialzones)

**Reads spatial transcriptomics input, computes tumor regions based on celltype annotations for tumor microenvironment analysis, and visualizes region assignments.**

## Features

- Computes **inside / interface / outside tumor regions** using spatial transcriptomics data (Xenium) 
- Utilizes **nearest neighbors** to define spatial context  
- Visualizes region assignments and regional compositions of cell types  
- Supports plotting gene expression across regions

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/murti-abhishek/spatialzones.git
cd spatialzones
pip install -e .
