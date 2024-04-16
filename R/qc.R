library(dplyr)
library(gridExtra)
library(ggplot2)
#' Generate QC Report Table for Seurat Object
#'
#' This function generates a QC report for a given Seurat object, allowing users to specify features and
#' output to a PDF file. It provides various statistical analyses including entropy, sparsity,
#' transcript counts per area, transcript counts per cell, and signal ratios. It also adjusts its behavior
#' based on the specified platform within the Seurat object (Xenium, CosMx, or Merscope).
#'
#' @param seu_obj A Seurat object, which must have properties like `path`, `platform`, and specific assay data.
#'                Must inherit from "Seurat".
#' @param features Optional; A character vector of features to include in the analysis.
#'                 Defaults to all features in `seu_obj` if `NULL`.
#' @param pdfFile A character string specifying the path and name of the output PDF file.
#'                Defaults to "QCReportTable.pdf".
#'
#' @return A data frame with the QC metrics for the provided Seurat object, which is also printed to the console.
#'         Additionally, a PDF file is generated containing a plot and table of the metrics.
#'
#' @details The function calculates various metrics such as entropy and sparsity using the `BioQC` and `coop` packages.
#'          It adjusts its processing logic based on the 'platform' attribute of the Seurat object to handle data from
#'          different technologies like Xenium, CosMx, or Merscope. It stops with an error if the object is not a Seurat object.
#'          The results are plotted using `ggplot2` and `gridExtra` for layout adjustments.
#'
#' @examples
#' # Assume 'seu_obj' is a preloaded Seurat object
#' genereateQCreport_table(seu_obj)
#'
#' @import dplyr
#' @import ggplot2
#' @import gridExtra
#' @importFrom BioQC entropy
#' @importFrom coop sparsity
#' @export
genereateQCreport_table <-
  function(seu_obj, features = NULL, pdfFile = "QCReportTable.pdf") {
    # function body as you provided
  }


genereateQCreport_table <-
  function(seu_obj, features = NULL, pdfFile = "QCReportTable.pdf") {
    if (!inherits(seu_obj, "Seurat")) {
      stop("Input must be a Seurat object.")
    }

    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    #conflicts_prefer("filter", "dplyr")
    path <- unique(seu_obj$path)
    platform <- unique(seu_obj$platform)
    tx_df <- readTxMeta(path, platform)
    # Extract metrics
    entropy_value = BioQC::entropy(as.matrix(seu_obj@assays$RNA$counts))
    sparsity_value = coop::sparsity(as.matrix(seu_obj@assays$RNA$counts))
    ncell = ncol(seu_obj)

    # More calculations
    tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
    neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)
    ratio = log10(tx_means) - log10(mean(neg_probe_means))
    mean_signal_ratio = mean(ratio)

    tx_count = colSums(seu_obj[["RNA"]]$counts[features,])
    mean_tx_norm = mean(tx_count / seu_obj$cell_area)
    tx_perarea = mean_tx_norm

    mean_tx = mean(colSums(seu_obj[["RNA"]]$counts[features,]))
    tx_percell = mean_tx

    tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
    neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)

    max_ratio <- log10(max(tx_means)) - log10(mean(neg_probe_means))


    if(platform == "Xenium"){
      #tx_df <- filter(tx_df, features %in% feature_name)
      tx_df <- tx_df[tx_df$feature_name %in% features, ]
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

    #max_vals <- matrixStats::rowMaxs(as.matrix(seu_obj[["RNA"]]$counts[features,]))
    #max_detection <-  max_vals
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      ncell = ncell,
      entropy_value = round(entropy_value, digits = 3),
      sparsity_value = round(sparsity_value, digits = 3),
      tx_perarea =  tx_perarea,
      tx_percell = tx_percell,
      cell_tx_fraction = cell_tx_fraction ,
      mean_ratio = mean_signal_ratio,
      max_ratio = max_ratio
      #max_detection= max_detection,

    )

    Metric_table = data.frame(
      Metric = c("Entropy", "Sparsity", "Tx per Area","Cell tx Fraction","Mean Singal Ratio","Max Ratio"),
      Value = c(entropy_value, sparsity_value, tx_perarea,cell_tx_fraction,mean_signal_ratio,max_ratio )
    )
    Metric_table$Metric <- factor(Metric_table$Metric, levels = Metric_table$Metric)

    pdf(pdfFile, width = 12, height = 4)

    p <- ggplot( Metric_table, aes(x = Metric, y = Value)) +
      geom_point(size = 4) +
      theme_minimal() +
      labs(title = "Single Sample Metrics", x = "", y = "Value")
    library(gridExtra)
    res_table <- tableGrob(res, rows = NULL)
    Metric_table <- tableGrob(Metric_table)
    grid.arrange(p,  Metric_table, ncol = 2, heights = c(5/6, 1/6))
    library(grid)
    grid.newpage()
    res_table <- tableGrob(res, rows = NULL)
    grid.draw(res_table)
    dev.off()

    print(res)
  }


genereateQCreport_table(seu_obj)
