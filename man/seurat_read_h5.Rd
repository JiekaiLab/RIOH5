% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratIO.R
\name{seurat_read_h5}
\alias{seurat_read_h5}
\title{H5 to Seuart object}
\usage{
seurat_read_h5(file = NULL, assay.name = NULL)
}
\arguments{
\item{file}{The h5 file}

\item{assay.name}{'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.}
}
\description{
Read h5 and converted h5 to the seurat object
}
