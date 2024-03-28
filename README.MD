# SpatialTouchstone

Welcome to `SpatialTouchstone`, an innovative tool designed to facilitate quality control for spatialomics technology data.

## Overview

Spatialomics technologies are revolutionizing how we understand biological tissues, enabling the mapping of molecular compositions at unprecedented spatial resolutions. However, with great power comes great complexity. The interpretation of spatial omics data requires robust quality control measures to ensure the reliability and reproducibility of the findings. `SpatialTouchstone` is developed to meet this need, providing a comprehensive suite of quality control metrics specifically tailored for spatial omics data.

## Functions - Internal use

### Functions in utils.R

#### readSpatial

' @title readSpatial. ' @description ' readSaptial reads in data from either Xenium, CosMx, of MERSCOPE and outputs a seurat object with some common metadata e for downstream comparison. ' @details ' Regardless of platform, data is stored in an assay named "RNA" for convenient and for each platform, this table will be used by subsequent functions. ' @param sample_id A unique identifier for the sample being read. ' @param path The file path to the directory containing the data files. ' @param platform A character string indicating the platform of the data, Must be one of "Xenium", "CosMx", or "MERSCOPE". ' @return A Seurat object containing the spatial data, with additional metadata and embeddings as applicable per platform. ' '@importFrom Seurat CreateSeuratObject' @importFrom data.table fread ' @export

#### readTxMeta

' @title readTxMeta.
' This function reads transcriptome metadata from a specified path, depending on the platform specified.
' Currently supports 'Xenium' and 'CosMx' platforms. For 'Xenium', it reads 'transcripts.csv.gz' and renames
' the 'feature_name' column to 'target' for consistency. For 'CosMx', it reads the appropriate transcriptome
' file matched by '*tx_file.csv.gz'. 'Merscope' platform support is under development.
' @description
'  It reads in the transcript localization/metadata table.
' @details
' For each platform, This table will be used by subsequent functions.
' @param path The file path from which to read the transcriptome metadata.
' @param platform  The platform for which the transcriptome metadata is being read. Valid options are
'        'Xenium', 'CosMx', and 'Merscope'. Note: 'Merscope' is currently not supported.
' @return Returns a data table with the transcriptome metadata. For 'Xenium', the 'feature_name' column is
'         renamed to 'target'. For 'CosMx', the transcriptome file matching '*tx_file.csv.gz' is read.
'         No return value for 'Merscope' as it stops execution with an error message.
' @export
' @importFrom data.table fread setnames

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

### Plot function in utils.R

-   plotSampleLabel
-   plotPanelSize
-   plotCellCount
-   plotTxPerCell
-   plotTxPerArea
-   plotTxPerNuc
-   plotTxPerCellNorm
-   plotFractionTxInCell
-   plotSignalRatio
-   plotMeanExpression
-   plotMaxExpression
-   plotMECR
-   plotMorans
-   plotSilhouette
-   plotSparsity
-   plotEntropy

## plotMetrics

## GetMetrics
