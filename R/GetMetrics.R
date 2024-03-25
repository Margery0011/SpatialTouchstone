#' @title getMetrics.
#' @description This function will take a path of where each folder is a sample.
#' @details
#'  This function processes a collection of spatial assay samples located in a specified path, where each folder corresponds to a separate sample.
#' It computes a range of spatial metrics, including mean expression, transcripts per cell, transcripts per area, signal-to-noise ratio, and more.
#' The function is highly flexible, allowing the user to specify which metrics to calculate.
#' Additionally, it can return the processed data in a Seurat object format if needed.
#' @param path  path to the directory containing sample folders. Each folder should be a separate sample.
#' @param mecr Logical; Default is TRUE, calculates the mean expression contrast ratio (MECR).
#' @param autocorr Logical; Default is TRUE, calculates Moran's I autocorrelation statistic.
#' @param clusterSilhouette Logical; Default is TRUE, calculates cluster silhouette widths.
#' @param sparsity Logical; Default is TRUE, calculates sparsity of the data.
#' @param entropy Logical; Default is TRUE, calculates Shannon entropy.
#' @param returnSeurat Logical; Default is FALSE, if TRUE, it returns a list with the original Seurat objects along with the calculated metrics.
#' @param sample_meta An optional data frame specifying sample metadata. If not provided, the function attempts to infer this information from the directory structure.
#' @return Depending on the value of `returnSeurat`, this function returns either a list containing both the original Seurat objects and the metrics data frame, or just the metrics data frame.
getMetrics <- function(path = NULL, mecr = TRUE, autocorr = TRUE,
                       clusterSilhouette = TRUE, sparsity = TRUE,
                       entropy = TRUE, returnSeurat = FALSE,
                       sample_meta = NULL) {
    # source("utils.R")
    # Create a data frame with information about samples found in path
    if(is.null(sample_meta)) {
        sample_meta <- data.frame(
        sample_id = list.files(path),
        path = paste0(path, list.files(path)),
        platform = stringr::word(list.files(path), 1, 1, sep="_")
        )
    }
    # Load objects from table
    tryCatch({
    obj_list <- list()
    for(i in 1:nrow(sample_meta)){
        obj_list[[i]] <- readSpatial(sample_id = sample_meta[i, "sample_id"],
        path = sample_meta[i, "path"],
        platform = sample_meta[i, "platform"])
        print(paste0("Done loading sample ", sample_meta$sample_id[i], ' ... ') )
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "")})
        print('Done loading all samples... ')
    ## Create list with probe set size and total cells - stored in a list
       metrics_list <- list(
        sample_meta = sample_meta,
    # Probe set size
      panel_size = data.frame(
      sample_id = sample_meta$sample_id,
      platform = sample_meta$platform,
      value = unlist(lapply(obj_list, nrow))
      ) ,
    # Total cells
    cell_count = data.frame(
    sample_id = sample_meta$sample_id,
    platform = sample_meta$platform,
    value = unlist(lapply(obj_list, ncol))
    ))
    ## apply functions in utils - keep adding to the metrics_list object
    metrics_list[['tx_per_cell']] <- do.call(rbind, lapply(obj_list, getTxPerCell))
    print("Done calculating transcripts per cell ... ")

    metrics_list[['tx_per_um2']] <- do.call(rbind, lapply(obj_list, getTxPerArea))
    print("Done calculating transcripts per um2 ... ")

    metrics_list[['tx_fraction_in_cell']] <- do.call(rbind, lapply(obj_list, getCellTxFraction))
    print("Done calculating transcripts per um2 ... ")

    # Transcripts per nucleus
    metrics_list[['tx_per_nuc']] <- do.call(rbind, lapply(obj_list, getTxPerNuc))
    print("Done calculating transcripts per nucleus ... ")

    # Transcript per cell (normalized by probe set size)
    metrics_list[['tx_per_cell_norm']] <- metrics_list[['tx_per_cell']]

    metrics_list[['signal_ratio']] <- do.call(rbind, lapply(obj_list, getMeanSignalRatio))
    print("Done calculating Mean signal-noise ratio ... ")


    # Expression
    metrics_list[['mean_expression']] <- do.call(rbind, lapply(obj_list, getMeanExpression))
    print("Done calculating Mean expression ... ")
    # MECR
    if(mecr == TRUE) {
        metrics_list[['mecr']] <- do.call(rbind, lapply(obj_list, getMECR))
        print("Done calculating MECR ... ")
    } else {print('Skipping MECR ... ')}

    # Moran's I
    if(autocorr == TRUE) {
        metrics_list[['morans']] <- do.call(rbind, lapply(obj_list, getMorans))
        print("Done calculating Moran's I autocorrelation ... ")
    } else {print("Skipping Moran's I autocorrelation ...  ")}

    # Silhouette Width
    if(clusterSilhouette == TRUE) {
        metrics_list[['silhouette']] <- do.call(rbind, lapply(obj_list, getSilhouetteWidth))
        print("Done calculating cluster silhouette width ... ")
    } else {print("Skipping cluster silhouette width ...  ")}

    #Sparsity
    if(sparsity == TRUE) {
        metrics_list[['sparsity']] <- do.call(rbind, lapply(obj_list, getSparsity))
        print("Done calculating Sparsity ... ")
    } else {print("Skipping sparsity ...  ")}


    #Entropy
    if(entropy == TRUE) {
        metrics_list[['entropy']] <- do.call(rbind, lapply(obj_list, getEntropy))
        print("Done calculating Shannon-Entropy ... ")
    } else {print("Skipping Shannon-Entropy ...  ")}

    if(returnSeurat == TRUE) {
        return(list(obj_list, metrics_list))
    } else {
        return(metrics_list)
  }
}
