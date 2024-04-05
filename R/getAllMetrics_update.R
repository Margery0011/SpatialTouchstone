# This function takes a data.frame and computes all metrics from files
# The data frame should have at least: one column called sample_id, platform, expMat, tx_file, and cell_meta.
# platform should be either "Xenium" or "CosMx", expMat = path to expression file, if platform is Xenium, this path should be to cell_feature_matrix folder,
# tx_file = path to tx file, cell_meta = path to metadata

source("R/utils_update.R")

getAllMetrics <- function(df_samples) {

  gStartTime <- Sys.time()

  tryCatch({
    start.time <- Sys.time()
    print("Calculating number of cells ")
    df_samples$nCells <- NA
    df_samples$nCells <- mapply(getNcells, expMat = df_samples$expMat, platform = df_samples$platform)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Specificity
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Specificity ")
    df_samples$specificityFDR <- NA
    df_samples$specificityFDR <- mapply(getGlobalFDR, tx_file = df_samples$tx_file, platform = df_samples$platform, cellSegMeta =  df_samples$cell_meta)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Tx per Cell
  tryCatch({
    start.time <- Sys.time()
    print("Calculating number of transcripts per cell ")
    df_samples$TxPerCell <- NA
    df_samples$TxPerCell <- mapply(getTxPerCell, expMat = df_samples$expMat, platform = df_samples$platform)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Tx per Area
  tryCatch({
    start.time <- Sys.time()
    print("Calculating number transcripts per area of segmented cell ")
    df_samples$TxPerArea <- NA
    df_samples$TxPerArea <- mapply(getTxPerArea, platform = df_samples$platform, cellSegMeta =  df_samples$cell_meta)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  # Tx per Nucleus
  tryCatch({
    start.time <- Sys.time()
    print("Calculating number transcripts per nuclei ")
    df_samples$TxPerNuc <- NA
    df_samples$TxPerNuc <- mapply(getTxPerNuc, tx_file = df_samples$tx_file, platform = df_samples$platform)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Signal to noise ratio
  tryCatch({
    start.time <- Sys.time()
    print("Calculating signal to noise ratio ")
    df_samples$SigNoiseRatio <- NA
    df_samples$SigNoiseRatio <- mapply(getMeanSignalRatio, expMat = df_samples$expMat, platform = df_samples$platform)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Fraction of txs in cells
  tryCatch({
    start.time <- Sys.time()
    print("Calculating fraction of transcripts in segmented cells ")
    df_samples$CellTxFraction <- NA
    df_samples$CellTxFraction <- mapply(getCellTxFraction, tx_file = df_samples$tx_file, platform = df_samples$platform)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})



  # MECR
  tryCatch({
    start.time <- Sys.time()
    print("Calculating MECR (Mutually Exclusive Co-Expression Rate) ")
    df_samples$MECR <- NA
    df_samples$MECR <- mapply(getMECR, expMat = df_samples$expMat, platform = df_samples$platform)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  # Sparsity
  tryCatch({
    start.time <- Sys.time()
    print("Calculating sparsity ")
    df_samples$sparsity <- NA
    df_samples$sparsity <- mapply(getSparsity, expMat = df_samples$expMat, platform = df_samples$platform)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Entropy
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Shannon entropy ")
    df_samples$entropy <- NA
    df_samples$entropy <- mapply(getEntropy, expMat = df_samples$expMat, platform = df_samples$platform)
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  gEndTime <- Sys.time()
  print(paste0("total time: ", round(gEndTime - gStartTime , digits = 2) ))
  return(df_samples)

}
