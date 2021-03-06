# --- scRNAw-seq analysis object H5 IO

# library
library(Seurat)
library(hdf5r)
library(Matrix)
library(SingleCellExperiment)
# library(monocle3)

#--- read h5 file

#' H5 to scRNAs-seq analysis object
#'
#' Read h5 and converted h5 to the scRNA-seq analysis object
#' @param target.object Denotes which object to load. Defualt is "seurat". Available options are:
#' \itemize{
#'   \item "seurat": converted the h5 file to "seurat object".
#'   \item "singlecellexperiment": converted the h5 file to "singlecellexperiment object".
#'   \item "monocle": converted the h5 file to "monocle3 object".
#' }
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type. Defualt is "RNA". Available options are:
#' \itemize{
#'   \item "RNA": this is the scRNA-seq data.
#'   \item "spatial": this is the spatial data.}
#' @return The single cell analysis mainstream software
#' @export
#'
read_h5 <- function(file, target.object = 'seurat', assay.name = 'RNA'){
  if(target.object == 'seurat'){
    data <- seurat_read_h5(file = file, assay.name = assay.name)
  }
  else if(target.object == 'singlecellexperiment'){
    data <- sce_read_h5(file = file, assay.name = assay.name)
  }
  else if(target.object == 'monocle'){
    data <- cds_read_h5(file = file, assay.name = assay.name)
  }
  else{
    stop('The output object must be specified')
  }
  return(data)
}


#--- write h5

#' The scRNAs-seq analysis objec to H5
#'
#' Write h5 and  the scRNA-seq analysis object converted to h5
#' @param object.type Denotes which object to save.Available options are:
#' \itemize{
#'   \item "seurat": converted the "seurat object" to the h5 file.
#'   \item "singlecellexperiment": converted the "singlecellexperiment object" to the h5 file.
#'   \item "monocle": converted the "monocle3 object" to the h5 file.
#' }
#' @param data The scRNA-seq analysis object data.
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type. Defualt is "RNA". Available options are:
#' \itemize{
#'   \item "RNA": this is the scRNA-seq data.
#'   \item "spatial": this is the spatial data.}
#' @param save.graphs Default is False , determing whether to save the graph(cell-cell similarity network).
#'   seurat graph is different from scanpy graph. Their relationship are set {"distances": "knn", "connectivities": "snn"} roughly.
#' @export
#'
write_h5 <- function(data, file, object.type = 'seurat', assay.name = 'RNA', save.graphs = FALSE){
  if(object.type == 'seurat'){
    seurat_write_h5(seurat = data, file = file, assay.name = assay.name, save.graphs = save.graphs)
  }
  if(object.type == 'singlecellexperiment'){
    sce_write_h5(sce = data, file = file, assay.name = assay.name)
  }
  if(object.type == 'monocle'){
    cds_write_h5(cds = data, file = file, assay.name = assay.name)
  }
}



