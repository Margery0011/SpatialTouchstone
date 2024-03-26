# SpatialTouchstone

Welcome to `SpatialTouchstone`, an innovative tool designed to facilitate quality control for spatialomics technology data.

## Overview

Spatialomics technologies are revolutionizing how we understand biological tissues, enabling the mapping of molecular compositions at unprecedented spatial resolutions. However, with great power comes great complexity. The interpretation of spatial omics data requires robust quality control measures to ensure the reliability and reproducibility of the findings. `SpatialTouchstone` is developed to meet this need, providing a comprehensive suite of quality control metrics specifically tailored for spatial omics data.

## Functions - Internal use

### Functions in utils.R

#### readSpatial

#### readTxMeta

It reads in the transcript localization/metadata table.

#### getGlobalFDR

Calculate Specificity as Global False Discovery Rate (FDR)

#### getTxPerCell

Calculates the average number of transcripts per cell for the given features.

#### getTxPerArea

#### getTxPerNuc

#### getMeanExpression

#### getMeanSignalRatio

#### getCellTxFraction

#### getMaxRatio

#### getMaxDetection

#### getMECR

#### getMorans

#### getSilhouetteWidth

#### getSparsity

#### getEntropy

###  Plot function in utils.R

- plotSampleLabel
- plotPanelSize 
- plotCellCount 
- plotTxPerCell
- plotTxPerArea
- plotTxPerNuc
- plotTxPerCellNorm 
- plotFractionTxInCell 
- plotSignalRatio 
- plotMeanExpression
- plotMaxExpression
- plotMECR 
- plotMorans 
- plotSilhouette 
- plotSparsity 
- plotEntropy
