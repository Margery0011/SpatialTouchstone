## This function take a list with the results from the function getMetrics to plot the metrics gathered
#' @title plotMetrics.
#' @param metrics_list A list containing the metrics results.
#' @param PlotAutocorr Logical; if `TRUE`, plots the autocorrelation metric. Default is `TRUE`.
#' @param PlotSparsity Logical; if `TRUE`, plots the sparsity metric;Default is `TRUE`.
#' @param PlotEntropy Logical; if `TRUE`, plots the entropy metric; Default is `TRUE`.
#' @param PlotClusterSilhouette Logical; if `TRUE`, plots the cluster silhouette metric; Default is `TRUE`.
#' @param ncol Integer; the number of columns in the plot grid; Default is 14.
#' @param nrow Integer; the number of rows in the plot grid; Default is 1.
#' @param rel_widths Numeric vector; relative widths of columns in the plot grid.
#' @return A plot object created with cowplot::plot_grid.
#' @importFrom cowplot plot_grid
#' @export

metrics_list <- list(
  sample_meta = sample_meta,
  getCellTxFraction_test_1 = getCellTxFraction_test_1,
  getTxPerNuc_test_1 = getTxPerNuc_test_1,
  etMeanSignalRatio_test_1 = getMeanSignalRatio_test_1,
  getMeanExpression_test_1 = getMeanExpression_test_1,
  getMorans_test = getMorans_test,
  etSparsity_test_1 = getSparsity_test_1,
  getEntropy_test_1 = getEntropy_test_1,
  getSilhouetteWidth_test = getSilhouetteWidth_test,
  tx_per_cell = getCellTxFraction_test_1, # Assuming this is the right mapping
  tx_per_nuc = getTxPerNuc_test_1,
  tx_per_cell_norm = getTxPerNuc_test_1, # Placeholder, adjust according to actual data
  tx_fraction_in_cell = getCellTxFraction_test_1, # Assuming repeated, using first instance
  signal_ratio = getMeanSignalRatio_test_1,
  mean_expression = getMeanExpression_test_1,
  sparsity = getSparsity_test_1,
  entropy = getEntropy_test_1,
  silhouette = getSilhouetteWidth_test,
  morans = getMorans_test # Assuming Moran's I test for spatial autocorrelation
)


plotMetrics <- function(metrics_list = NULL, PlotAutocorr = T, PlotSparsity = T,
                        PlotEntropy = T,
                        PlotClusterSilhouette = T,
                        ncol = 14, nrow = 1,
                        rel_widths = c(0.75, 0.4, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,0.3,0.75, 0.75, 0.8, 0.3, 0.8, 0.8)) {
# PLOT
    p0 <- plotSampleLabel(metrics_list[['sample_meta']])
    #p1 <- plotPanelSize(metrics_list[["panel_size"]])
    #p2 <- plotCellCount(metrics_list[["cell_count"]])
    p3 <- plotTxPerCell(metrics_list[["tx_per_cell"]])
    #p4 <- plotTxPerArea(metrics_list[["tx_per_um2"]])
    p5 <- plotTxPerNuc(metrics_list[["tx_per_nuc"]])
    p6 <- plotTxPerCellNorm(metrics_list[["tx_per_cell_norm"]])
    p7 <- plotFractionTxInCell(metrics_list[["tx_fraction_in_cell"]])
    p8 <- plotSignalRatio(metrics_list[["signal_ratio"]])
    p9 <- plotMeanExpression(metrics_list[["mean_expression"]])
    #p10 <- plotMECR(metrics_list[["mecr"]])
  if(PlotAutocorr == TRUE) {
    p11 <- plotMorans(metrics_list[["morans"]])
    }
  if(PlotSparsity == TRUE) {
    p13 <- plotSparsity(metrics_list[["sparsity"]])
    }
  if(PlotEntropy == TRUE) {
    p14 <- plotEntropy(metrics_list[["entropy"]])
    }
  if(PlotClusterSilhouette == TRUE) {
    p12 <- plotSilhouette(metrics_list[["silhouette"]])
  }
    p <- cowplot::plot_grid(p0, #p1, #p2,
                            p3, #p4,
                            p5, p6,  p7, p8, p9,
                            #p10,
                            p11, p12,p13, p14,
        ncol=ncol, align='h',
        rel_widths = rel_widths , nrow = nrow)

  return(p)
}
