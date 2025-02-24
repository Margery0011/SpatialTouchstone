% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_update_final.R
\name{getMaxDetection}
\alias{getMaxDetection}
\title{getMaxDetection}
\usage{
getMaxDetection(
  seu_obj = NULL,
  features = NULL,
  expMat = "path_to_expMat",
  platform = NULL
)
}
\arguments{
\item{seu_obj}{A Seurat object.}

\item{features}{An optional vector of gene identifiers for which to analyze the distribution of maximal
expression values. If NULL, the calculation encompasses all genes in the dataset.}

\item{expMat}{Path to exprMatrix file.}

\item{platform}{The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
and 'Merscope'.}
}
\value{
A data frame summarizing the maximal expression values across the specified features or the
entire dataset, which can be used to analyze the upper limits of detection and expression within the
sample.

A data frame.
}
\description{
Calculate the distribution of maximal values.
}
\details{
This function identifies the maximal expression values across the specified set of features
(or all features if none are specified) within the dataset, illustrating the upper bounds of gene
expression. Such information is crucial for assessing the dataset's dynamic range and the sensitivity
of detection methods used in the experiment. The function aggregates these maximal values and presents
them in a data frame, facilitating further analysis of the expression distribution and detection
efficiency across different genes.
}
