library(Seurat)
library(dplyr)
library(ggplot2)
library(shadowtext)
library(scales)
library(cowplot)
library(data.table)
library(Matrix)
library(matrixStats)
library(SingleCellExperiment)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(bluster)
library(BiocParallel)
library(scales)

#######
# I/O
#######

# Functions assume data is stored in specific formats:
# Xenium: Typical Xenium bundle; CosMx: Flat CSV (TAP style); MERSCOPE: Not implemented yet

#####
# readSpatial() reads in data from either Xenium, CosMx, of MERSCOPE.
# It outputs a seurat object with some common metadata e for downstream comparison.
# Regardless of platform, data is stored in an assay named "RNA" for convenient
# function
# If Seurat = F, then store sample in a list, where each item of the list is a data.table (one for Tx information, one for expression, etc)
readSpatial <- function(sample_id, path, platform=NULL, seurat=FALSE){
  print(paste0("Reading: ", sample_id))

  # If platform not specified, try to guess which tech it is (from folder name)
  if(is.null(platform)) {
    xenium <- grep('*[X-x]enium*|*10x*', sample_id)
    cosmx <- grep('*[C-c]os[M-m]x*|*[N-n]anostring*', sample_id)
    mersc <- grep('*[M-m]erscope*|*MERSCOPE*|*[V-v]izgen*',sample_id)

    if(length(xenium) == 0 & length(cosmx) == 0 & length(mersc) == 0) {
      platform = 'Unknown'
    }

    if(length(xenium) > 0) { platform = 'Xenium'}
    if(length(cosmx) > 0) { platform = 'CosMx'}
    if(length(mersc) > 0) { platform = 'Merscope'}
  }

  ## Read and store tables as list outside of seurat
  if(seurat==FALSE) {
    # create empty list
    obj_list <- list()

    if(platform == 'Xenium') {
      obj_list[[sample_id]] <- list()
      obj_list[[sample_id]][['expMatrix']] <- Matrix::readMM(file.path(path, 'cell_feature_matrix/matrix.mtx.gz'))
      cols <- data.table::fread(file.path(path, 'cell_feature_matrix/barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(path, 'cell_feature_matrix/features.tsv.gz'), header = F)
      rownames(obj_list[[sample_id]][['expMatrix']]) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(obj_list[[sample_id]][['expMatrix']]) <- cols$V1 ## this is the barcodes of cells

      obj_list[[sample_id]][['TxMatrix']] <- data.table::fread(file.path(path, 'transcripts.csv.gz'))
    }

    if(platform == 'CosMx') {
      obj_list[[sample_id]] <- list()
      pattern <- "exprMat_file.csv.gz$|exprMat*"
      file_name <- list.files(path = path, pattern = pattern, full.names = TRUE)
      obj_list[[sample_id]][['expMatrix']] <- data.table::fread(file_name)
      pattern <- "*tx_file.csv.gz$|*tx.csv$|tx_file_unique.csv.gz$"
      file_name <- list.files(path = path, pattern = pattern, full.names = TRUE)
      obj_list[[sample_id]][['TxMatrix']] <- data.table::fread(file_name[1])
    }

    # if(platform == 'Merscope') {
    #   obj_list[[sample_id]] <- list()
    #   pattern <- "exprMat_file.csv.gz$"
    #   file_name <- list.files(path = path, pattern = pattern, full.names = TRUE)
    #   obj_list[[sample_id]][['expMatrix']] <- data.table::fread(file_name)
    #   obj_list[[sample_id]][['TxMatrix']] <- data.table::fread(file.path(path, 'transcripts.csv.gz'))
    # }
    #

    return(obj_list)
  }


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
    seu_obj[["tissue"]] <- CreateDimReducObject(m , key="Tissue_", assay="RNA")


  } else if(platform == "Merscope"){
    print("Working on support!")
    stop()

  } else{
    print("Not a supported platform")
    stop()

  }

  return(seu_obj)
}

#####
# readTxMeta() simply reads in the transcript localization/metadata table
# for each platform. This table will be used by subsequent functions
readTxMeta <- function(path, platform, sample_meta=NULL){
  if(platform == "Xenium"){
    df <- data.table::fread(file.path(path, "transcripts.csv.gz"))
    ## change feature_name to target - to keep consistency ##
    setnames(df, "feature_name", "target")

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

## number of cells
getNcells <- function(seu_obj = NULL, expMat = 'path_to_expMat', platform = NULL) {
  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      ncell <- ncol(Matrix::readMM(file.path(expMat, 'matrix.mtx.gz')))
    }
    if(platform == 'CosMx') {
      ncell <- nrow(data.table::fread(expMat))
    }
  }

  if(!is.null(seu_obj)) {
    ncell <- ncol(seu_obj)
  }

  return(ncell)

}


### Specificity as Global FDR ###
getGlobalFDR <- function(seu_obj = NULL, features = NULL, tx_file ='path_to_txFile',platform = NULL) {

  if(is.null(seu_obj)) {
    tx_df <- data.table::fread(tx_file)

    if(platform == 'Xenium') {
      setnames(tx_df, "feature_name", "target") ## changing the colname to target to keep it consistent for all techs
    }

    if(platform == 'CosMx') {
      tx_df <- data.table::fread(tx_file)

    }

    if(platform == 'Merscope') {

    }
  }


  if(!is.null(seu_obj)) {
    ## Obj seurat should have a column in @meta.data with path to Tx file
    path <- unique(seu_obj$path)
    # same for Platform
    platform <- unique(seu_obj$platform)

    # Read Tx localization data
    tx_df <- readTxMeta(path, platform)
  }


  ## Get probes that are Negative control or 'blank' barcodes ##
  negProbes <- tx_df$target[grep('Neg*|SystemControl*|Blank*|BLANK*', tx_df$target)]
  allGenes <- unique(tx_df[!target %in% negProbes, target]) ## list of unique genes (non control or blank probes) in panel

  # create table with expression per each gene in panel (adding N of Txs for each gene in object allGenes) - p.s This will take the expression of Txs outside of assigned cells
  expTableAll  <- tx_df[, .(Count = .N), by = target]
  expTable <- expTableAll[expTableAll$target %in% allGenes, ]
  expNeg <- sum(expTableAll[!expTableAll$target %in% allGenes, ]$Count) ## sum of all Negative control or blank or unassigned barcodes (i.e non specific)

  numGenes <- length(expTable$target)
  numNeg <- length(expTableAll[!expTableAll$target %in% allGenes, ]$target)

  expTable$FDR <- 1
  for(i in allGenes) {
    fdr = (expNeg / (expTable[expTable$target %in% i, ]$Count + expNeg) ) * (numGenes / numNeg) * 1/100
    expTable[target == i, FDR := fdr]
  }

  if(is.null(features)) {
    return(mean(expTable$FDR))
  } else {
    return(mean(expTable$FDR[expTable$target %in% features]))
  }

}




### Transcripts per cell
#features can be explicitly defined. Defaults to all targets
## expMat = specify where the "cell_feature_matrix" folder is IF Xenium. IF cosmx - path to exp matrix
getTxPerCell <- function(seu_obj = NULL, features=NULL, expMat = 'path_to_exprMatrix',
                         platform){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      exp <- Matrix::readMM(file.path(expMat, 'matrix.mtx.gz'))
      cols <- data.table::fread(file.path(expMat, 'barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(expMat, 'features.tsv.gz'), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells

    }

    if(platform == 'Merscope') {

    }

    if(is.null(features)) {
      features <- rownames(exp)
      # remove non specific probes
      features <- features[-grep('Unassigned*|NegControl*|BLANK*|SystemControl*', features)]
    }

    # Calculate average N of Txs per cell
    mean_tx <- mean(colSums(exp[features, ]))

    return(mean_tx)

  }

  if(!is.null(seu_obj)) {
    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }


    mean_tx <- mean(colSums(seu_obj[["RNA"]]$counts[features,]))
    return(mean_tx)
  }

}




### Transcripts per um2
## If features are specified, have to specify where the Tx file is
getTxPerArea <- function(seu_obj = NULL, features=NULL,
                         platform, cellSegMeta = 'path_to_cellMeta', tx_file = NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      cell_meta <- data.table::fread(file.path(cellSegMeta))

      mean_tx_norm <- mean(cell_meta$total_counts / cell_meta$cell_area)

      # If features are specified - have to calculate differently
      if(!is.null(features)) {
        ## Have to read the transcripts file
        tx_df <- data.table::fread(tx_file)
        # subset the tx file with only existing assigned cells and features
        tx_df <- tx_df[tx_df$cell_id %in% cell_meta$cell_id,]
        tx_df <- tx_df[tx_df$feature_name %in% features,]
        # count number of features per cell - have to merge cell_meta and tx_df - to add the area information
        tx_df <- merge(tx_df, cell_meta, by = 'cell_id')
        tx_df <- tx_df %>% group_by(cell_id,cell_area) %>% tally()

        mean_tx_norm <- mean(tx_df$n / tx_df$cell_area)
      }

    }

    if(platform == 'CosMx') {
      cell_meta <- data.table::fread(file.path(cellSegMeta))
      mean_tx_norm <- mean(cell_meta$nCount_RNA / cell_meta$Area.um2)

      if(!is.null(features)) {
        ## Have to read the transcripts file
        tx_df <- data.table::fread(tx_file)
        # subset the tx file with only existing assigned cells and features
        tx_df <- tx_df[cell_ID != 0 & target %in% features]
        # count number of features per cell - have to merge cell_meta and tx_df - to add the area information
        tx_df <- merge(tx_df, cell_meta, by = 'cell')
        tx_df <- tx_df %>% group_by(cell,Area.um2) %>% tally()

        mean_tx_norm <- mean(tx_df$n / tx_df$Area.um2)
      }

    }

    if(platform == 'Merscope') {

    }

    return(mean_tx_norm)

  }

  if(!is.null(seu_obj)) {
    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    tx_count <- colSums(seu_obj[["RNA"]]$counts[features,])
    mean_tx_norm <- mean(tx_count / seu_obj$cell_area)


    return(mean_tx_norm)
  }


}

### Transcripts per nucleus
getTxPerNuc <- function(seu_obj=NULL, features=NULL, tx_file = NULL, platform = NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      tx_df <- data.table::fread(tx_file)

      if(is.null(features)) {
        # subset the tx file with only existing assigned cells and features
        # remove neg control probes
        negProbes <- unique(tx_df$feature_name[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', tx_df$feature_name)])
        # number of txs in nucleus
        nTx_nuc <- dim(tx_df[cell_id != 'UNASSIGNED' & overlaps_nucleus == 1 & !feature_name %in% negProbes])[1]
        # number of cells
        nCells <- length(unique(tx_df$cell_id))
      }
      if(!is.null(features)) {
        nTx_nuc <- dim(tx_df[cell_id != 'UNASSIGNED' & overlaps_nucleus == 1 & !feature_name %in% negProbes & feature_name %in% features])[1]

      }


    }

    if(platform == 'CosMx') {
      tx_df <- data.table::fread(tx_file)

      if(is.null(features)) {
        # subset the tx file with only existing assigned cells and features
        # remove neg control probes
        negProbes <- unique(tx_df$target[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', tx_df$target)])
        # number of txs in nucleus
        nTx_nuc <- nrow(tx_df[cell_ID != 0 & CellComp == 'Nuclear' & !target %in% negProbes])
        # number of cells
        nCells <- length(unique(tx_df$cell[tx_df$cell_ID != 0]))
      }
      if(!is.null(features)) {
        nTx_nuc <- nrow(tx_df[cell_ID != 0 & CellComp == 'Nuclear' & !target %in% negProbes & target %in% features])

      }

    }

    if(platform == 'Merscope') {

    }

    return(nTx_nuc / nCells)

  }

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
    tx_df <- filter(tx_df, cell_id %in% colnames(seu_obj) &
                      overlaps_nucleus == 1 &
                      features %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())


  } else if(platform == "CosMx"){
    tx_df$cell_id <- paste(tx_df$cell_ID, tx_df$fov, sep="_")
    tx_df <- tx_df %>%
      filter(cell_id %in% colnames(seu_obj) &
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
# This function will take in a exp matrix or a seurat object, and return a df with average expression per probe,
# You can specify which genes, or default, all genes.
# For xenium, when specifying the exp matrix - should be the path to the cell_feature folder
getMeanExpression <- function(seu_obj = NULL, features=NULL, expMat = 'path_to_expMatrix',
                              platform=NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      exp <- Matrix::readMM(file.path(expMat, 'matrix.mtx.gz'))
      cols <- data.table::fread(file.path(expMat, 'barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(expMat, 'features.tsv.gz'), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells
      # remove neg control
      exp <- exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)),]
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      ## remove first 2 columns - usually FOV and Cell_ID information
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells

    }

    if(platform == 'Merscope') {

    }

    avg_exp_df <- as.data.frame(rowMeans(exp))
    colnames(avg_exp_df) <- 'MeanExpression'

    if(is.null(features)) {
      return(avg_exp_df)
    } else {
      return(subset(avg_exp_df, rownames(avg_exp_df) %in% features))
    }

  }



  if(!is.null(seu_obj)) {
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

}

### log-ratio of mean gene counts to mean neg probe counts (signal to noise ratio)
getMeanSignalRatio <- function(seu_obj=NULL, features=NULL, platform=NULL, expMat = 'path_to_expMat'){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      exp <- Matrix::readMM(file.path(expMat, 'matrix.mtx.gz'))
      cols <- data.table::fread(file.path(expMat, 'barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(expMat, 'features.tsv.gz'), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells

    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      ## remove first 2 columns - usually FOV and Cell_ID information
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }

    noise <- exp[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)), ]
    exp <- exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)),]
    #ratio <- mean( log10(rowMeans(exp + .1)) - log10(rowMeans(noise + .1)) )

    if(is.null(features)) {
      return( suppressWarnings(mean( log10(rowMeans(exp + .1)) - log10(rowMeans(noise + .1)) )))
    } else {
      return(suppressWarnings(mean( log10(rowMeans(exp[rownames(exp) %in% features,] + .1)) - log10(rowMeans(noise + .1)) )))
    }


  }

  if(!is.null(seu_obj)) {
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


}

### Fraction of transcripts in cells
# This function takes in a seurat object or a Tx file and return the fraction of txs in assigned cells compared to total
getCellTxFraction <- function(seu_obj=NULL, features=NULL, tx_file = 'path_to_tx',
                              platform = NULL,path){

  if(is.null(seu_obj)) {

    tx_df <- data.table::fread(tx_file)
    total_tx_count <- nrow(tx_df)

    if(platform == 'Xenium') {
      if(is.null(features)) {
        unassigned_tx_count <- sum(tx_df$cell_id == 'UNASSIGNED')
      } else {
        tx_df <- tx_df[tx_df$feature_name %in% features,]
        total_tx_count <- nrow(tx_df)
        unassigned_tx_count <- sum(tx_df$cell_id == 'UNASSIGNED')

      }

    }

    if(platform == 'CosMx') {
      if(is.null(features)) {
        unassigned_tx_count <- sum(tx_df$CellComp == '')

      } else {
        tx_df <- tx_df[tx_df$target %in% features,]
        total_tx_count <- nrow(tx_df)
        unassigned_tx_count <- sum(tx_df$CellComp == '')

      }

    }

    return( (total_tx_count - unassigned_tx_count) / total_tx_count )

  }


  if(!is.null(seu_obj)) {
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
      tx_df <- tx_df[tx_df$target%in% features, ]
      total_tx_count <- nrow(tx_df)
      unassigned_tx_count <- sum(tx_df$cell_id == "UNASSIGNED")

      cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count

    } else if(platform == "CosMx"){
      tx_df <- filter(tx_df, target %in% features)
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



}

##### Dynamic Range
# Log-ratio of highest mean exp vs. mean noise
getMaxRatio <- function(seu_obj = NULL, features=NULL, expMat ='path_to_expMat', platform = NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      exp <- Matrix::readMM(file.path(expMat, 'matrix.mtx.gz'))
      cols <- data.table::fread(file.path(expMat, 'barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(expMat, 'features.tsv.gz'), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      ## remove first 2 columns - usually FOV and Cell_ID information
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }

    tx_means <- rowMeans(exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)),])
    neg_probe_means <- rowMeans(exp[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)), ])

    if(is.null(features)) {
      return( log10(max(tx_means)) - log10(mean(neg_probe_means)) )

    } else {
      return( log10(max(tx_means[features])) - log10(mean(neg_probe_means)) )
    }

  }


  if(!is.null(seu_obj)) {

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
}

# Distribution of maximal values
getMaxDetection <- function(seu_obj = NULL, features=NULL, expMat ='path_to_expMat', platform = NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      exp <- Matrix::readMM(file.path(expMat, 'matrix.mtx.gz'))
      cols <- data.table::fread(file.path(expMat, 'barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(expMat, 'features.tsv.gz'), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      ## remove first 2 columns - usually FOV and Cell_ID information
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }

    max_vals <- data.frame(matrixStats::rowMaxs(as.matrix(exp)))
    colnames(max_vals) <- 'MaxValue'

    if(is.null(features)) {
      return(max_vals)
    } else {
      return(subset(max_vals, rownames(max_vals) %in% features))
    }

  }



  if(!is.null(seu_obj)) {
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

}

##### Mutually Exclusive Co-expression Rate (MECR) Implementation
getMECR <- function(seu_obj=NULL, expMat = 'path_to_expMat', platform = NULL) {
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

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      exp <- Matrix::readMM(file.path(expMat, 'matrix.mtx.gz'))
      cols <- data.table::fread(file.path(expMat, 'barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(expMat, 'features.tsv.gz'), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      ## remove first 2 columns - usually FOV and Cell_ID information
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }

    genes <- intersect(rownames(exp), rownames(marker_df))
    mtx <- as.matrix(exp[genes,])

  }

  if(!is.null(seu_obj)) {

    genes <- intersect(rownames(seu_obj), rownames(marker_df))
    mtx <- as.matrix(seu_obj[['RNA']]$counts[genes, ])


  }

  coexp.rates <- c()
  #print(paste0("Marker count: ", length(genes)))
  if (length(genes) > 25) { genes <- sample(genes, 25) }
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


  return(round(mean(coexp.rates), digits=3))

}

##### Distribution spatial autocorrelation
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
  sce <- SingleCellExperiment(list(counts=seu_obj[["RNA"]]$counts[features,]),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)

  colGraph(sfe, "knn20") <- findSpatialNeighbors(sfe, method = "knearneigh",
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
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)

  #Nearest neighbor
  colGraph(sfe, "knn20") <- findSpatialNeighbors(sfe, method = "knearneigh",
                                                 dist_type = "idw", k = 20,
                                                 style = "W")
  #Moran's I
  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = MulticoreParam(8))

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
getSparsity <- function(seu_obj=NULL, features = NULL, expMat = 'path_to_expMat', platform = NULL) {

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      exp <- Matrix::readMM(file.path(expMat, 'matrix.mtx.gz'))
      cols <- data.table::fread(file.path(expMat, 'barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(expMat, 'features.tsv.gz'), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }

    if(is.null(features)) {
      return(coop::sparsity(as.matrix(exp)))
    } else {
      return(coop::sparsity(as.matrix(exp[rownames(exp) %in% features, ])))
    }
  }

  if(!is.null(seu_obj)) {
    value = coop::sparsity(as.matrix(seu_obj@assays$RNA$counts))
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=round(value, digits=3)
    )
    return(res)
  }

}


getEntropy <- function(seu_obj=NULL, features = NULL, expMat = 'path_to_expMat', platform = NULL) {

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      exp <- Matrix::readMM(file.path(expMat, 'matrix.mtx.gz'))
      cols <- data.table::fread(file.path(expMat, 'barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(expMat, 'features.tsv.gz'), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      exp <- exp[, -c(1:2)]
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }

    if(is.null(features)) {
      return(BioQC::entropy(as.matrix(exp)))
    } else {
      return(BioQC::entropy(as.matrix(exp[rownames(exp) %in% features, ])))
    }
  }


  if(!is.null(seu_obj)) {
    value = BioQC::entropy(as.matrix(seu_obj@assays$RNA$counts))
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=round(value, digits=3)
    )
    return(res)
  }

}


#######
# Plotting
#######
# All plots assume input is a tidy data frame with the following columns:
# 1) sample_id
# 2) platform
# 3) value (based on what is being plotted--from functions above)

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

# plotPanelSize <- function(df){
#   p <- ggplot(df, aes(x="", y=sample_id)) +
#     geom_point(shape=21, color='black', alpha=0.8, stroke=1,
#                aes(size=value, fill=value)) +
#     geom_shadowtext(color = "black", size = 4, #fontface = "bold",
#                     bg.colour = "white", bg.r = .2,
#                     aes(label=scales::comma(value))) +
#     scale_fill_gradientn(colours=viridis::mako(100)) +
#     xlab("Panel size") + ylab("") +
#     scale_size(range = c(6,12)) +
#     scale_x_discrete(position='top',
#                      labels = c("")) +
#     theme_classic() +
#     theme(
#       legend.position='none',
#       axis.text.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.title.x = element_text(size=12),
#       axis.line.y = element_blank(),
#       axis.ticks.y = element_blank()
#     )
#
#   return(p)

#}

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

#plotTxPerCellIntersect <-

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

#plotTxPerCell_Intersect <- function()

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
                  labels = label_log(digits = 2)) +
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



