# SpatialTouchstone

Welcome to `SpatialTouchstone`, an innovative tool designed to facilitate quality control for spatialomics technology data.

## Overview

Spatialomics technologies are revolutionizing how we understand biological tissues, enabling the mapping of molecular compositions at unprecedented spatial resolutions. However, with great power comes great complexity. The interpretation of spatial omics data requires robust quality control measures to ensure the reliability and reproducibility of the findings. `SpatialTouchstone` is developed to meet this need, providing a comprehensive suite of quality control metrics specifically tailored for spatial omics data.
# readSpatial() reads in data from either Xenium, CosMx, of MERSCOPE.
# It outputs a seurat object with some common metadata e for downstream comparison.
# Regardless of platform, data is stored in an assay named "RNA" for convenient
# function
