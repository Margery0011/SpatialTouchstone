

#* @apiTitle SpatialMultiQC
#* @apiDescription Quality control for probe-based Spatial Assays


#* Calculate number of transcripts per cell
#* @param exp_table Expression table - each column as a cell and rows as features - csv format
#* @param features Features of interest - if NULL will use every feature available
#* @post /TxPerCell

function(exp_table, features=NULL){
  if(is.null(features)){
    features <- rownames(exp_table)
  } else{
    features <- features
  }

  exp_table <- read.csv(exp_table, header = T, row.names = 1)

  mean_tx <- mean(colSums(exp_table))

  return(mean_tx)
}



### Transcripts per um2


#* Calculate Transcripts per um2
#* @param exp_table Expression table - each column as a cell and rows as features - csv format
#* @param features Features of interest - if NULL will use every feature available
#* @post /TxPerArea
function(seu_obj,
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

#* Calculate Transcripts per nucleus
#* @param exp_table Expression table - each column as a cell and rows as features - csv format
#* @param features Features of interest - if NULL will use every feature available
#* @post /TxPerNuc
function(seu_obj,
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

#
### Per Probe Mean Expression

#* Calculate Mean expression
#* @param exp_table Expression table - each column as a cell and rows as features - csv format
#* @param features Features of interest - if NULL will use every feature available
#* @post /MeanExpression
#*
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
#* Calculate log-ratio of mean gene counts to mean neg probe counts
#* @param exp_table Expression table - each column as a cell and rows as features - csv format
#* @param features Features of interest - if NULL will use every feature available
#* @post /MeanSignalRatio

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
