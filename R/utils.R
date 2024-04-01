
#######
# I/O
#######



globalVariables(c("target", "FDR", "cell_id","overlaps_nucleus","CellComp","platform","value","sample_id","type","samples","column","squish"))
globalVariables(c(":=", "."))



#' Reads and preprocesses spatial transcriptomics data from a specified path for a given platform.
#' The function supports 'Xenium' and 'CosMx' platforms, performing platform-specific loading and
#' preprocessing steps. This includes loading the data, annotating it with sample metadata, and
#' processing cell metadata. For 'Xenium', it adds additional cell metadata and tissue coordinates
#' as an embedding for custom plotting. For 'CosMx', it fixes assays to separate targeting and
#' non-targeting probes, adds additional cell metadata, and includes tissue coordinates as an embedding.
#' 'Merscope' platform support is under development.
#' @title readSpatial.
#' @description
#' readSaptial reads in data from either Xenium, CosMx, of MERSCOPE and outputs a seurat object with some common metadata e for downstream comparison.
#' @details
#' Regardless of platform, data is stored in an assay named "RNA" for convenient and for each platform, this table will be used by subsequent functions.
#' @param sample_id Identifier for the sample being processed.
#' @param path The file path from which to load the spatial transcriptomics data.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
##' @return Returns a Seurat object containing the loaded and processed spatial transcriptomics data,
#'         including sample metadata, cell metadata, and embeddings for tissue coordinates. The
#'         structure and content of the returned object vary depending on the platform.
#' @importFrom Seurat CreateSeuratObject
#' @importFrom data.table fread
#' @examples
#' seu_obj <- readSpatial(sample_id = "sample123", path = "path/to/xenium/data", platform = "Xenium")
#' seu_obj <- readSpatial(sample_id = "sample456", path = "path/to/cosmx/data", platform = "CosMx")
#' @export
#' @keywords readSpatialdata
readSpatial <- function(sample_id, path, platform){
  print(paste0("Reading: ", sample_id))
  if(platform == "Xenium"){
    print("Loading Xenium data")
    seu_obj <- LoadXenium(path, assay = "RNA")
    seu_obj@meta.data$sample_id <- sample_id
    seu_obj@meta.data$platform <- platform
    seu_obj@meta.data$path <- path #used in some functions to pull tx data

    ##### Cell Metadata
    print("Getting additional cell metadata")
    cell_meta <- data.table::fread(
      file.path(path, "cells.csv.gz")
    )
    #Set up a few defined metadata columns
    seu_obj@meta.data$cell_area <- cell_meta$cell_area
    seu_obj@meta.data$nucleus_area <- cell_meta$nucleus_area
    seu_obj@meta.data$transcript_counts <- cell_meta$transcript_counts
    seu_obj@meta.data$negprobe_counts <- cell_meta$control_probe_counts

    ### Add tissue coordinates as an embedding for custom plotting
    print("Adding tissue coordinates as embedding")
    coords <- GetTissueCoordinates(seu_obj)
    coords <- as.matrix(coords[,1:2])
    colnames(coords) <- c("Tissue_1", "Tissue_2")
    rownames(coords) <- colnames(seu_obj)
    seu_obj[["tissue"]] <- CreateDimReducObject(coords, key="Tissue_", assay="RNA")

  } else if(platform == "CosMx"){
    print("Loading CosMx data")
    seu_obj <- LoadNanostring(path, fov="fov")
    seu_obj$sample_id <- sample_id
    seu_obj@meta.data$platform <- platform
    seu_obj@meta.data$path <- path #used in some functions to pull tx data

    #Fix assays to separate targeting and non-targeting probes
    ## Add negative control probe assay
    sys_probes <- grep("SystemControl", rownames(seu_obj), value=T)
    neg_probes <- grep("Negative", rownames(seu_obj), value=T)
    seu_obj[["ControlProbe"]] <- CreateAssayObject(
      counts = seu_obj[["Nanostring"]]$counts[neg_probes,]
    )
    ## Make "Nanostring" assay
    tx_probes <- rownames(seu_obj)[!rownames(seu_obj) %in% c(sys_probes, neg_probes)]
    seu_obj[["RNA"]] <- CreateAssayObject(
      counts = seu_obj[["Nanostring"]]$counts[tx_probes,]
    )
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj[["Nanostring"]] <- NULL

    ##### Cell Metadata
    print("Getting additional cell metadata")
    cell_meta <- data.table::fread(
      file.path(path, list.files(path, pattern="*metadata_file.csv.gz"))
    )
    #It's excessive, but we'll add all metadata into the object
    seu_obj@meta.data <- cbind(seu_obj@meta.data, cell_meta)
    seu_obj@meta.data$fov <- factor(paste0("FOV", seu_obj@meta.data$fov))
    seu_obj@meta.data$cell_area <- seu_obj$Area.um2
    seu_obj@meta.data$transcript_counts <- seu_obj$nCount_RNA
    seu_obj@meta.data$negprobe_counts <- seu_obj$nCount_ControlProbe

    ### Add tissue coordinates as an embedding for custom plotting
    print("Adding tissue coordinates as embedding")
    coords <- data.frame(
      Tissue_1 = cell_meta$CenterY_global_px,
      Tissue_2 = cell_meta$CenterX_global_px
    )
    coords <- as.matrix(coords[,1:2])
    colnames(coords) <- c("Tissue_1", "Tissue_2")
    rownames(coords) <- colnames(seu_obj)
    seu_obj[["tissue"]] <- CreateDimReducObject(coords, key="Tissue_", assay="RNA")


  } else if(platform == "Merscope"){
    print("Working on support!")
    stop()

  } else{
    print("Not a supported platform")
    stop()

  }

  return(seu_obj)
}


#' This function reads transcriptome metadata from a specified path, depending on the platform specified.
#' Currently supports 'Xenium' and 'CosMx' platforms. For 'Xenium', it reads 'transcripts.csv.gz' and renames
#' the 'feature_name' column to 'target' for consistency. For 'CosMx', it reads the appropriate transcriptome
#' file matched by '*tx_file.csv.gz'. 'Merscope' platform support is under development.
#' @title readTxMeta.
#' @description
#'  It reads in the transcript localization/metadata table.
#' @details
#' For each platform, This table will be used by subsequent functions.
#' @param path The file path from which to read the transcriptome metadata.
#' @param platform  The platform for which the transcriptome metadata is being read. Valid options are
#'        'Xenium', 'CosMx', and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @return Returns a data table with the transcriptome metadata. For 'Xenium', the 'feature_name' column is
#'         renamed to 'target'. For 'CosMx', the transcriptome file matching '*tx_file.csv.gz' is read.
#'         No return value for 'Merscope' as it stops execution with an error message.
#' @export
#' @importFrom data.table fread setnames

readTxMeta <- function(path, platform){
  if(platform == "Xenium"){
    df <- data.table::fread(file.path(path, "transcripts.csv.gz"))
    ## change feature_name to target - to keep consistency ##
    data.table::setnames(df, "feature_name", "target")

  } else if(platform == "CosMx"){
    df <- data.table::fread(file.path(path,
                                   list.files(path, pattern = "*tx_file.csv.gz")))
  } else if(platform == "Merscope"){
    print("Working on support!")
    stop()
  } else{
    print("Platform not supported")
    stop()
  }
}

#######
# QC
#######
#' @details
#' Computes the global FDR for specified features (or all features by default) in a Seurat object.
#' It leverages transcript localization data to calculate FDR based on the proportion of negative control
#' or blank barcode expressions compared to the expression of actual genes.
#' @title getGlobalFDR
#' @description
#' Calculate Specificity as Global False Discovery Rate (FDR).
#' @param seu_obj A Seurat object containing spatial data, including a path and platform attribute.
#' @param features An optional vector of feature names (gene names) for which to calculate the FDR.
#'        If NULL (default), FDR is calculated for all features in the seu_obj.
#' @export
#' @importFrom data.table fread .N
#' @importFrom Seurat CreateSeuratObject
#' @return A data frame with columns for sample_id, platform, and the
#' calculated mean FDR across the specified features.
#' @examples
#'
#' res <- getGlobalFDR(seu_obj)
#' res <- getGlobalFDR(seu_obj, features = c("Gene1", "Gene2"))
getGlobalFDR <- function(seu_obj, features=NULL) {
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)

  # Read Tx localization data
  tx_df <- readTxMeta(path, platform)

  ## Get probes that are Negative control or 'blank' barcodes ##
  negProbes <- tx_df$target[grep('Neg*|SystemControl*|Blank*|BLANK*', tx_df$target)]
  allGenes <- unique(tx_df[!target %in% negProbes, target]) ## list of unique genes (non control or blank probes) in panel

  # create table with expression per each gene in panel (adding N of Txs for each gene in object allGenes) - p.s This will take the expression of Txs outside of assigned cells
  expTableAll  <- tx_df[, .(Count = .N), by = target]
  expTable <- expTableAll[expTableAll$target %in% allGenes, ]
  expNeg <- sum(expTableAll[!expTableAll$target %in% allGenes, ]$Count) ## sum of all Negative control or blank or unassigned barcodes

  numGenes <- length(expTable$target)
  numNeg <- length(expTableAll[!expTableAll$target %in% allGenes, ]$target)

  expTable$FDR <- 1
  for(i in allGenes) {
    fdr = (expNeg / (expTable[expTable$target %in% i, ]$Count + expNeg) ) * (numGenes / numNeg) * 1/100
    expTable[target == i, FDR := fdr]
  }

  res <- data.frame(
    #sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value= mean(expTable$FDR)
  )
  return(res)

}



#' @title getTxPerCell.
#' @description Calculates the average number of transcripts per cell for the given features.
#' @details
#' This function calculates the mean number of transcripts per cell for a specified set of features (genes)
#' in a Seurat object. If no features are specified, the function defaults to using all targets available
#' within the RNA assay of the provided Seurat object. It's a useful metric for assessing the overall
#' transcriptional activity within the sampled cells.
#' @param seu_obj A Seurat object with RNA assays. The object must have 'sample_id' and 'platform' metadata attributes for identification
#'        and reporting purposes.
#' @param features An optional vector of feature names (e.g., gene symbols) to include in the calculation.
#'        Defaults to NULL, in which case the calculation uses all available features in the RNA assay
#'        of the Seurat object.
#' @return A data frame with columns 'sample_id', 'platform', and 'value', where 'value' represents the
#'         mean number of transcripts per cell calculated across the specified features. This output can
#'         be useful for comparative analysis across samples or experimental conditions.
#' @export
#' @import Seurat
getTxPerCell <- function(seu_obj, #features can be explicitly defined. Defaults to all targets
                         features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  mean_tx <- mean(colSums(seu_obj[["RNA"]]$counts[features,]))

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean_tx
  )
  return(res)
}

#' This function calculates the mean number of transcripts per unit area for specified features (genes) in a Seurat object.
#' If no features are specified, the function defaults to using all targets within the RNA assay of the provided Seurat object.
#' This calculation provides insight into the density of transcriptional activity relative to cell area, offering a normalized
#' measure of gene expression that accounts for differences in cell size.
#' @title getTxPerArea.
#' @description
#' Calculate Mean Transcripts per Unit Area
#' @param seu_obj A Seurat object.
#'        The object must have 'sample_id' and 'platform' metadata for identification and reporting. It is expected that
#'        `seu_obj$cell_area` contains the area information for each cell.
#' @param features An optional vector of feature names (e.g., gene symbols) for which to calculate transcripts per unit area.
#'        Defaults to NULL, which means the calculation uses all available features in the RNA assay of the Seurat object.
#' @return Returns a data frame with three columns: `sample_id`, `platform`, and `value`. The `value` column contains
#'         the computed mean number of transcripts per unit area for the selected features across all cells in the dataset.
#'         This data frame provides a concise summary of transcriptional activity normalized by cell size, which can be
#'         critical for downstream analyses, especially in studies where cell morphology and size are variable.
#' @import Seurat
#' @export
getTxPerArea <- function(seu_obj,
                         features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  tx_count <- colSums(seu_obj[["RNA"]]$counts[features,])
  mean_tx_norm <- mean(tx_count / seu_obj$cell_area)

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean_tx_norm
  )
  return(res)
}

### Transcripts per nucleus
#' @title getTxPerNuc.
#' @description
#' It calculates ranscripts per nucleus.
#' @details
#' This function calculates the mean number of transcripts localized within the nucleus
#' for a specified set of features (genes) across all cells in a given dataset. This measurement
#' is specific to spatial transcriptomics data where the distinction between nuclear and cytoplasmic
#' transcript localization can be made. The function is designed to work with different platforms by
#' reading transcript localization data using a helper function `readTxMeta`. If no features are
#' specified, the function defaults to using all available genes.
#' @param seu_obj A Seurat object containing spatial transcriptomics data with metadata attributes
#'        including 'path' and 'platform' to identify the location and type of data. The object is
#'        expected to be annotated with sample identifiers and must have been preprocessed to include
#'        information on cellular localization of transcripts.
#' @param features An optional vector of feature names (genes) for which to calculate the mean transcripts
#'        per nucleus. If NULL (default), the function calculates this metric for all genes.
#' @return A data frame with the columns: `sample_id`, `platform`, and `value`, where `value` represents
#'         the mean number of transcripts per nucleus across the specified features for the dataset.
#'         This provides a focused view of nuclear gene expression.
#'
#' @examples
#' `seu_obj` is a Seurat object prepared with the necessary metadata
#' res <- getTxPerNuc(seu_obj)
#' # For specific genes:
#' res  <- getTxPerNuc(seu_obj, features = c("GAPDH", "ACTB"))
#' @export
#' @import Seurat
#' @importFrom dplyr filter group_by summarize

getTxPerNuc <- function(seu_obj,
                        features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)

  # Read Tx localization data
  tx_df <- readTxMeta(path, platform)

  if(platform == "Xenium"){
    tx_df <- dplyr::filter(tx_df, cell_id %in% colnames(seu_obj) &
                      overlaps_nucleus == 1 &
                      features %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())


  } else if(platform == "CosMx"){
    tx_df$cell_id <- paste(tx_df$cell_ID, tx_df$fov, sep="_")
    tx_df <- tx_df %>%
      dplyr::filter(cell_id %in% colnames(seu_obj) &
               CellComp == "Nuclear" &
               target %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())

  } else if(platform == "Merscope"){
    print("Working on support")

  } else{
    print("Platform not supported")
  }

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean(tx_df$nuc_counts)
  )

  return(res)
}


### Per Probe Mean Expression
#' @title getMeanExpression.
#' @description
#' It calculated mean expression per probe.
#' @details
#' Computes the mean expression levels for a specified list of features (genes)
#' and control probes within a Seurat object. This function separately calculates the mean
#' expression for both target genes and control probes, then combines these results into a
#' single data frame. If no features are specified, it defaults to all genes in the RNA assay.
#' @param seu_obj A Seurat object.
#' This object must have an RNA assay for target genes and a ControlProbe assay for control probes.
#' @param features An optional vector of gene names for which to calculate mean expression.
#' If NULL, mean expression is calculated for all genes in the RNA assay.
#' @return A data frame with columns `target` (gene or control probe name), `value` (mean expression level),
#' `type` (indicating whether the row represents a Gene or Control), `platform`, and `sample_id`.
#' @export
#' @import Seurat
getMeanExpression <- function(seu_obj,
                              features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  target_df <- data.frame(
    target = features,
    value = rowMeans(seu_obj[["RNA"]]$counts[features,]),
    type = "Gene"
  )

  control_df <- data.frame(
    target = rownames(seu_obj[["ControlProbe"]]$counts),
    value = rowMeans(seu_obj[["ControlProbe"]]$counts),
    type = "Control"
  )

  res <- rbind(target_df, control_df)
  res$platform <- unique(seu_obj$platform)
  res$sample_id <- unique(seu_obj$sample_id)

  return(res)
}

### log-ratio of mean gene counts to mean neg probe counts
#' @title getMeanSignalRatio.
#' @description
#' It calculates log-ratio of mean gene counts to mean neg probe counts.
#' @details
#' Computes the log-ratio of the mean expression levels of specified genes (or all genes if none are specified)
#' to the mean expression levels of negative control probes within a Seurat object. This metric can provide insights into
#' the overall signal strength relative to background noise levels.
#' @param seu_obj A Seurat object, expected to have both 'RNA' and 'ControlProbe'
#' assays for calculating mean expression levels of genes and negative control probes, respectively.
#' @param features An optional vector of gene names for which to calculate the mean log-ratio. If NULL, the calculation
#' is performed for all genes present in the 'RNA' assay of the seu_obj.
#' @return A data frame with three columns: `sample_id`, `platform`, and `value`, where `value` represents the calculated
#' mean log-ratio of gene expression to negative control probe expression for the selected genes. This summary can be
#' used to assess the signal-to-noise ratio in the dataset.
#' @export

getMeanSignalRatio <- function(seu_obj,
                               features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
  neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)

  ratio <- log10(tx_means) - log10(mean(neg_probe_means))
  ratio <- mean(ratio)

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=ratio
  )

  return(res)
}

### Fraction of transcripts in cells
#' @title getCellTxFraction.
#' @description
#' It calculates fraction of transcripts in cells.
#' @details
#' The function works by reading transcript metadata for the given sample and platform, then filtering for the specified features (if any).
#' For each platform, it identifies transcripts that are not assigned to any cell ('UNASSIGNED' for Xenium, 'None' for CosMx, etc.) and calculates the fraction of transcripts that are assigned to cells.
#' This metric is an important indicator of the quality of spatial transcriptomics data, reflecting the efficiency of transcript capture and cell assignment.
#' @param seu_obj A Seurat object containing spatial transcriptomics data, expected to include metadata fields for `path` and `platform` to facilitate reading transcript metadata with `readTxMeta`.
#' @param features An optional vector of gene identifiers for which to calculate the transcript fraction. If NULL, the calculation is performed using all available features in the dataset.
#' @return A data frame containing the calculated fraction of transcripts in cells, along with `sample_id` and `platform` for context. This output provides insight into the efficiency of transcript capture and assignment within the given dataset.
#' @return A data frame with probe type (Gene or Control), the mean expression values, and the associated sample and platform identifiers.
#' @export
#' @import Seurat
getCellTxFraction <- function(seu_obj,
                              features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)

  tx_df <- readTxMeta(path, platform)

  if(platform == "Xenium"){
    #tx_df <- filter(tx_df, features %in% feature_name)
    tx_df <- tx_df[tx_df$feature_name %in% features, ]
    total_tx_count <- nrow(tx_df)
    unassigned_tx_count <- sum(tx_df$cell_id == "UNASSIGNED")

    cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count

  } else if(platform == "CosMx"){
    tx_df <- dplyr::filter(tx_df, target %in% features)
    total_tx_count <- nrow(tx_df)
    unassigned_tx_count <- sum(tx_df$CellComp == "None")

    cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count

  } else if(platform == "Merscope"){
    print("Working on support")

  } else{
    print("Platform not supported")
  }

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=cell_tx_fraction
  )

  return(res)
}

##### Dynamic Range
# Log-ratio of highest mean exp vs. mean noise
#' @title getMaxRatio.
#' @description
#' It calculateds the Log-ratio of highest mean exp vs. mean noise.
#' @details
#' The function identifies the maximum mean expression value among the specified features
#' (or all features if none are specified) and calculates its log-ratio to the mean expression value
#' of negative control probes. This log-ratio reflects the dynamic range of the dataset, indicating
#' the spread between the highest signal and background noise levels, which is critical for assessing
#' data quality and sensitivity of detection in spatial transcriptomics.
#' @param seu_obj A Seurat object including assays for target genes ('RNA')
#' and control probes ('ControlProbe').
#' @param features An optional vector of gene identifiers for which to perform the calculation. If NULL, the calculation
#' includes all genes in the RNA assay.
#' @return A data frame with `sample_id`, `platform`, and the calculated `value` of the log-ratio.
#' @export
#' @import Seurat
getMaxRatio <- function(seu_obj,
                        features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
  neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)

  ratio <- log10(max(tx_means)) - log10(mean(neg_probe_means))

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=ratio
  )

  return(res)
}

# Distribution of maximal values
#' @title getMaxDetection.
#' @description
#' It calculates the distribution of maximal values.
#' @details This function identifies the maximal expression values across the specified set of features
#' (or all features if none are specified) within the dataset, illustrating the upper bounds of gene
#' expression. Such information is crucial for assessing the dataset's dynamic range and the sensitivity
#' of detection methods used in the experiment. The function aggregates these maximal values and presents
#' them in a data frame, facilitating further analysis of the expression distribution and detection
#' efficiency across different genes.
#'
#' @param seu_obj A Seurat object.
#' @param features An optional vector of gene identifiers for which to analyze the distribution of maximal
#' expression values. If NULL, the calculation encompasses all genes in the dataset.
#' @return A data frame summarizing the maximal expression values across the specified features or the
#' entire dataset, which can be used to analyze the upper limits of detection and expression within the
#' sample.
#' @return A data frame.
#' @export
#' @import Seurat

getMaxDetection <- function(seu_obj,
                            features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  max_vals <- matrixStats::rowMaxs(as.matrix(seu_obj[["RNA"]]$counts[features,]))

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=max_vals,
    gene = features
  )

  return(res)
}

##### Mutually Exclusive Co-expression Rate (MECR) Implementation
#' @title getMECR.
#' @description
#' It calculates  Mutually Exclusive Co-expression Rate (MECR) Implementation.
#' @details
#' Based on a predefined set of markers and their associated cell types, this function computes the rate
#' at which pairs of markers are expressed in mutually exclusive patterns within cells. The calculation is limited
#' to markers found in the dataset and can be adjusted to focus on a subset by specifying it. The result is a
#' single MECR value that quantifies the overall mutually exclusive co-expression pattern, rounded to three decimal places.
#' @param  seu_obj A seurat object.
#' @return A data frame containing the `sample_id`, `platform`, and the computed
#' MECR value. The MECR value is rounded to three decimal places and represents
#' the average mutually exclusive co-expression rate across the selected markers.
#' @export
#' @import Seurat


getMECR <- function(seu_obj) {
  #This function comes from Hartman & Satija, bioRxiv, 2024
  #We are using a custom marker table. The original publication bases it on
  #scRNA-seq from matched tissue.
  marker_df <- data.frame(
    gene = c("EPCAM", "KRT19", "KRT8",
             "CD3E", "CD3D", "CD8A", "NKG7",
             "MS4A1", "CD79A",
             "PECAM1", "CLDN5", "VWF",
             "C1QA", "C1QB", "CD14", "FCGR3A", "ITGAX", "ITGAM",
             "PDGFRA", "DPT", "COL1A1",
             "MYH11", "ACTG2"),
    cell_type = c("Epithelial", "Epithelial", "Epithelial",
                  "T", "T", "T", "T",
                  "B", "B",
                  "Endo", "Endo", "Endo",
                  "Macro", "Macro", "Macro", "Macro", "Macro", "Macro",
                  "Fibro", "Fibro", "Fibro",
                  "Muscle", "Muscle")
  )
  rownames(marker_df) <- marker_df$gene

  coexp.rates <- c()
  genes <- intersect(rownames(seu_obj), rownames(marker_df))
  print(paste0("Marker count: ", length(genes)))
  if (length(genes) > 25) { genes <- sample(genes, 25) }
  mtx <- as.matrix(seu_obj[['RNA']]$counts[genes, ])
  for (g1 in genes) {
    for (g2 in genes) {
      if ((g1 != g2) && (g1 > g2) && (marker_df[g1, "cell_type"] != marker_df[g2, "cell_type"])) {
        c1 <- mtx[g1, ]
        c2 <- mtx[g2, ]
        coexp.rates <- c(
          coexp.rates,
          sum(c1 > 0 & c2 > 0) / sum(c1 > 0 | c2 > 0)) # >0 too liberal of an expression threshold?
      }
    }
  }

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(mean(coexp.rates), digits=3)
  )

  return(res)
}

##### Distribution spatial autocorrelation

#' @title getMorans.
#' @param  seu_obj A seurat object.
#' @param features A character vector specifying the features (genes or probes) to include in the analysis. If `NULL` (the default), all features in `seu_obj` are used.
#' @return A data frame.
#' @details
#' The function first processes gene-targeting probes by converting the input data into a `SingleCellExperiment` object, then into a `SpatialFeatureExperiment` object;
#' Spatial coordinates are specified and used to further process the data, including normalization and nearest neighbor graph construction. Subsequently, Moran's I is calculated for each feature using the `Voyager` package. The same process is repeated for control probes. The results for both gene-targeting and control probes are combined into a single data frame.
#' The spatial autocorrelation analysis provides Moran's I values for each feature, indicating the degree of spatial clustering. This is performed separately for gene-targeting probes and control probes, allowing for a comparison of spatial patterns in gene expression and background noise.
#' @return A data frame with Moran's I values for each feature (gene or control probe), along with sample ID, platform, gene/probe name, and type (Gene or Control). Each row corresponds to a feature, and columns include the sample identifier, platform used, Moran's I value, gene/probe name, and the type indicating whether it's a Gene-targeting probe or a Control probe.
#' @export
#' @import Seurat
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importMethodsFrom SingleCellExperiment counts
#' @importFrom SpatialExperiment toSpatialExperiment
#' @importFrom SpatialFeatureExperiment toSpatialFeatureExperiment
#' @importFrom Voyager runMoransI
#' @importFrom scater logNormCounts
#' @importFrom BiocParallel MulticoreParam


getMorans <- function(seu_obj,
                      features=NULL){
  #Requires SingleCellExperiment, SpatialFeatureExperiment, Voyager, scater

  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  #First run for gene-targeting probes
  print("Getting Moran's I for gene-targeting probes")
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=seu_obj[["RNA"]]$counts[features,]),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- SpatialExperiment::toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- SpatialFeatureExperiment::toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  SpatialFeatureExperiment::rowData(sfe)$means <- rowMeans(counts(sfe))
  SpatialFeatureExperiment::rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)

  SpatialFeatureExperiment::colGraph(sfe, "knn20") <- SpatialFeatureExperiment::findSpatialNeighbors(sfe, method = "knearneigh",
                                                 dist_type = "idw", k = 20,
                                                 style = "W")

  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = MulticoreParam(8))

  spatial_cor <- as.data.frame(rowData(sfe))

  targeting <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=spatial_cor[,3], #morans I
    gene = rownames(spatial_cor),
    type = "Gene"
  )

  #Now run for control probes
  print("Getting Moran's I for non-targeting probes")
  sce <- SingleCellExperiment(list(counts=seu_obj[["ControlProbe"]]$counts),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- SpatialExperiment::toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(SingleCellExperiment::counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)

  #Nearest neighbor
  SpatialFeatureExperiment::colGraph(sfe, "knn20") <- SpatialFeatureExperiment::findSpatialNeighbors(sfe, method = "knearneigh",
                                                 dist_type = "idw", k = 20,
                                                 style = "W")
  #Moran's I
  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = BiocParallel::MulticoreParam(8))

  spatial_cor <- as.data.frame(rowData(sfe))

  control <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=spatial_cor[,3], #morans I
    gene = rownames(spatial_cor),
    type = "Control"
  )

  res <- rbind(targeting, control)

  return(res)
}

##### Cluster evaluation: silhouette width
#' @title getSilhouetteWidth.
#' @description
#' It calculates silhouette width as cluster evaluation.
#' @details The silhouette width calculation involves preprocessing steps including normalization and scaling of the data, PCA for dimensionality reduction, and clustering.
#' The silhouette width is then calculated for a downsampled subset of cells to ensure computational efficiency.
#' This metric helps in assessing the cohesion and separation of the identified clusters, with higher values indicating better defined clusters.
#' @param seu_obj A Seurat object.
#' @return A data frame containing the `sample_id`, `platform`,
#' and the average silhouette width across all clusters, rounded to three decimal places.
#' @export
#' @import Seurat
#' @importFrom bluster approxSilhouette


getSilhouetteWidth <- function(seu_obj){
  print("Clustering data")
  seu_obj <- seu_obj %>%
    NormalizeData() %>%
    ScaleData()
  VariableFeatures(seu_obj) <- rownames(seu_obj)
  seu_obj <- seu_obj %>%
    RunPCA(verbose=F) %>%
    FindNeighbors(dims=1:10) %>%
    FindClusters(resolution=0.2)

  #Downsample to 100 cells per cluster for silhouette calculation
  seu_obj <- subset(seu_obj, downsample = 10000)

  silhouette <- bluster::approxSilhouette(
    Embeddings(seu_obj, 'pca')[,1:10],
    clusters = seu_obj$seurat_clusters
  )

  silhouette <- as.data.frame(silhouette)

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(mean(silhouette$width), digits=3)
  )

}

##### Sparsity calculation #####
#Show the sparsity (as a count or proportion) of a matrix.
#For example, .99 sparsity means 99% of the values are zero. Similarly, a sparsity of 0 means the matrix is fully dense.
#' @title getSparsity.
#' @description
#' It shows the sparsity (as a count or proportion) of a matrix.
#' @details
#' For example, .99 sparsity means 99% of the values are zero. Similarly, a sparsity of 0 means the matrix is fully dense.
#' @param seu_obj A Seurat object.
#' @return A data frame.
#' @export
#' @importFrom BioQC entropy
#' @import Seurat
#' @importFrom coop sparsity
getSparsity <- function(seu_obj) {
  value = coop::sparsity(as.matrix(seu_obj@assays$RNA$counts))
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(value, digits=3)
  )
  return(res)
}

#' @title getEntropy.
#' @description Computes the entropy of the RNA count matrix within a Seurat object.
#' @details
#' Entropy is used to quantify the diversity or uniformity of gene expression across the dataset. This function calculates the entropy of the RNA count matrix, reflecting the distribution of gene expression levels.
#' A higher entropy value suggests a more uniform distribution across genes, while a lower value indicates concentration in a smaller number of genes.
#' This measure can provide insights into the complexity of the cellular composition and the heterogeneity within the sample.
#' @param seu_obj A Seurat object containing RNA count data.
#' @return A data frame with the `sample_id`, `platform`, and the entropy value of the RNA count matrix, rounded to three decimal places.
#' @export
#' @importFrom BioQC entropy
getEntropy <- function(seu_obj) {
  value = BioQC::entropy(as.matrix(seu_obj@assays$RNA$counts))
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(value, digits=3)
  )
  return(res)
}



#######
# Plotting
#######
# All plots assume input is a tidy data frame with the following columns:
# 1) sample_id
# 2) platform
# 3) value (based on what is being plotted--from functions above)

#' @title plotSampleLabel
#' @param sample_meta A data frame containing at least two columns: `sample_id` and `platform`.
#' @return A ggplot object representing the sample labels colored by platform.
#' @import ggplot2
#' @export

plotSampleLabel <- function(sample_meta){
  df <- data.frame(
    samples = sample_meta$sample_id,
    platform = sample_meta$platform
  )

  p <- ggplot(df, aes(x="", y=samples)) +
    geom_text(aes(label=samples, color=platform), size=5, hjust=0.5) +
    scale_color_manual(values = c("#59C134", "#14B3E6")) +
    theme_void() + theme(legend.position='none')
  return(p)
}

#' @title plotPanelSize
#' @param df A data frame expected to contain `sample_id` and `value`, where `value` represents the panel size.
#' @return A ggplot object visualizing the panel sizes of samples.
#' @export
#' @import ggplot2
#' @importFrom viridis mako
#' @importFrom scales comma
#' @importFrom shadowtext geom_shadowtext
plotPanelSize <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    shadowtext::geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("Panel size") + ylab("") +
    scale_size(range = c(6,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)

}

#' @title plotCellCount.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotCellCount <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_text(aes(label=scales::comma(value), color=platform), size=5, hjust=0.5) +
    scale_color_manual(values = c("#59C134", "#14B3E6")) +
    scale_x_discrete(position='top') +
    xlab("Cell count") + ylab("") +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  return(p)
}

#' @title plotTxPerCell.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotTxPerCell <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("Per cell") + ylab("") +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotTxPerArea.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotTxPerArea <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per cell\num^2") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}


#' @title plotTxPerNuc.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotTxPerNuc <- function(df){
  df$column <- ""
  p <- ggplot(df, aes(x=column, y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per\nnucleus") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotTxPerCellNorm.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotTxPerCellNorm <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per cell\nper target") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotFractionTxInCell.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotFractionTxInCell <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90', stroke=1) +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -.05) +
    xlab("Fraction Tx in Cells") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0),
                       breaks=c(0, 0.5, 1),
                       labels = c(0, 0.5, 1),
                       limits=c(0,1)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}


#' @title plotSignalRatio.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotSignalRatio <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90') +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -.05) +
    xlab("Mean log10-ratio\nexpression over noise") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}
#' @title plotMeanExpression.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom scales label_log
plotMeanExpression <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=type),
                position = position_jitterdodge(jitter.width=0.25)) +
    geom_boxplot(color="black",
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 aes(fill=type)) +
    scale_fill_manual(values=c("lightgrey", "firebrick")) +
    scale_colour_manual(values=c("grey20", "firebrick")) +
    xlab("Mean gene\n detection per cell") + ylab("") +
    scale_x_log10(position='top', expand = c(0,0),
                  labels = scales :: label_log(digits = 2)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotMaxExpression.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotMaxExpression <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=platform)) +
    geom_boxplot(color="black", fill="lightgrey",
                 alpha=0.5, outlier.size=0, outlier.colour = NA) +
    scale_colour_manual(values = c("#59C134", "#14B3E6")) +
    xlab("Maximal gene\ndetection per cell") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0),
                  limits = c(0,50), oob=squish) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

#' @title plotMECR.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotMECR <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9, "YlOrRd"),
                         limits=c(0.01)) +
    xlab("MECR") + ylab("") +
    scale_size(range = c(7,12), limits = c(0, 0.1)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotMorans.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotMorans <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=type),
                position = position_jitterdodge(jitter.width=0.25)) +
    geom_boxplot(color="black",
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 aes(fill=type)) +
    scale_fill_manual(values=c("lightgrey", "firebrick")) +
    scale_colour_manual(values=c("grey20", "firebrick")) +
    xlab("Mean gene\n spatial autocorrelation\n(Moran's I)") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotSilhouette.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotSilhouette <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90') +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -0.001) +
    xlab("Mean\nsilhouette width\n(Louvain res=0.2)") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotSparsity.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotSparsity <- function(df){
    p <- ggplot(df, aes(x="", y=sample_id)) +
          geom_point(shape=21, color='black', alpha=0.8, stroke=1,
          aes(size=value, fill=value)) +
          geom_shadowtext(color = "black", size = 4, #fontface = "bold",
          bg.colour = "white", bg.r = .2,
          aes(label=scales::comma(value))) +
          scale_fill_gradientn(colours=viridis::mako(100)) +
          xlab("Sparsity") + ylab("") +
          scale_size(range = c(7,12)) +
          scale_x_discrete(position='top',
          labels = c("")) +
          theme_classic() +
          theme(
              legend.position='none',
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size=10),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank()
             )
    return(p)
}

#' @title plotEntropy.
#' @param df A data frame.
#' @return A ggplot object.
#' @export
#' @import ggplot2

plotEntropy <- function(df){
    p <- ggplot(df, aes(x="", y=sample_id)) +
          geom_point(shape=21, color='black', alpha=0.8, stroke=1,
          aes(size=value, fill=value)) +
        geom_shadowtext(color = "black", size = 4, #fontface = "bold",
          bg.colour = "white", bg.r = .2,
          aes(label=scales::comma(value))) +
          scale_fill_gradientn(colours=viridis::mako(100)) +
          xlab("Entropy") + ylab("") +
          scale_size(range = c(7,12)) +
          scale_x_discrete(position='top',
          labels = c("")) +
          theme_classic() +
          theme(
            legend.position='none',
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size=10),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank()
          )
  return(p)
}



