# SWJGL 

## Analyzing Omics Data to Unveil Dependence Patterns Across Neural Layers

Code for exploring gene expression networks across cortical layers using spatially-resolved transcriptomics data from human dorsolateral prefrontal cortex.

**Authors:** Alex Cecchetto¹, Davide Forcina², Mariafrancesca Patalano¹  
¹Department of Statistical Sciences, University of Padova  
²Department of Industrial Engineering, University of Naples Federico II

### Results Overview

![Network analysis results for each cortical layer](images/layers.pdf)

Analyzes gene expression networks across six grey matter layers and white matter of the human brain using two spatially-aware methods:
- **SWJGL**: Spatially Weighted Joint Graphical Lasso
- **SWFLSA**: Spatially Weighted Fused Lasso Signal Approximation

Spatial weights are computed using the Hausdorff distance between layers.

### Usage

See `scripts/main.R`.

This repository is related to the Statistical Models course (PhD course - unipd).
