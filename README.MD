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


' @details
' Computes the global FDR for specified features (or all features by default) in a Seurat object.
' It leverages transcript localization data to calculate FDR based on the proportion of negative control
' or blank barcode expressions compared to the expression of actual genes.
' @title getGlobalFDR
' @description
' Calculate Specificity as Global False Discovery Rate (FDR).
' @param seu_obj A Seurat object containing spatial data, including a path and platform attribute.
' @param features An optional vector of feature names (gene names) for which to calculate the FDR.
'        If NULL (default), FDR is calculated for all features in the seu_obj.
' @export
' @importFrom data.table fread .N
' @importFrom Seurat CreateSeuratObject
' @return A data frame with columns for sample_id, platform, and the
' calculated mean FDR across the specified features.
' @examples

#### getTxPerCell

Calculates the average number of transcripts per cell for the given features.

' @title getTxPerCell.
' @description Calculates the average number of transcripts per cell for the given features.
' @details
' This function calculates the mean number of transcripts per cell for a specified set of features (genes)
' in a Seurat object. If no features are specified, the function defaults to using all targets available
' within the RNA assay of the provided Seurat object. It's a useful metric for assessing the overall
' transcriptional activity within the sampled cells.
' @param seu_obj A Seurat object containing single-cell RNA sequencing data, including predefined assays
'        like 'RNA'. The object must have 'sample_id' and 'platform' metadata attributes for identification
'        and reporting purposes.
' @param features An optional vector of feature names (e.g., gene symbols) to include in the calculation.
'        Defaults to NULL, in which case the calculation uses all available features in the RNA assay
'        of the Seurat object.
' @return A data frame with columns 'sample_id', 'platform', and 'value', where 'value' represents the
'         mean number of transcripts per cell calculated across the specified features. This output can
'         be useful for comparative analysis across samples or experimental conditions.
' @export
' @import Seurat

#### getTxPerArea

' This function calculates the mean number of transcripts per unit area for specified features (genes) in a Seurat object.
' If no features are specified, the function defaults to using all targets within the RNA assay of the provided Seurat object.
' This calculation provides insight into the density of transcriptional activity relative to cell area, offering a normalized
' measure of gene expression that accounts for differences in cell size.
' @title getTxPerArea.
' @description
' Calculate Mean Transcripts per Unit Area
' @param seu_obj A Seurat object containing single-cell RNA sequencing data, with an RNA assay and cell area information.
'        The object must have 'sample_id' and 'platform' metadata for identification and reporting. It is expected that
'        `seu_obj$cell_area` contains the area information for each cell.
' @param features An optional vector of feature names (e.g., gene symbols) for which to calculate transcripts per unit area.
'        Defaults to NULL, which means the calculation uses all available features in the RNA assay of the Seurat object.
' @return Returns a data frame with three columns: `sample_id`, `platform`, and `value`. The `value` column contains
'         the computed mean number of transcripts per unit area for the selected features across all cells in the dataset.
'         This data frame provides a concise summary of transcriptional activity normalized by cell size, which can be
'         critical for downstream analyses, especially in studies where cell morphology and size are variable.
' @import Seurat
' @export

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
