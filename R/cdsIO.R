# --- monocle(cds) H5 IO

# library
library(hdf5r)
library(Matrix)
library(monocle)

# --- monocle read the h5 file

#' H5 to monocle object
#'
#' Read h5 and converted h5 to the monocle object
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @return monocle object
#'
cds_read_h5 <- function(file=NULL, assay.name = NULL){
  if(!file.exists(file)){
    stop('No such file or directory')
  }
  h5 <- H5File$new(filename = file, mode = 'r')
  tryCatch({
    data <- h5_to_cds(h5 = h5, assay.name = assay.name)
  },
  error = function(e) {
    print(e)
  },
  finally = {
    h5$close_all()
  }
  )
  return(data)
}

#' H5 to the monocle3
#'
#' H5 file is converted to the monocle3 obejct
#' @param h5 The h5 file in R
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#'
h5_to_cds <- function(h5, assay.name='RNA'){
  if(h5attr(h5, 'assay_name') == assay.name){
    cds_list <- list()
    for(i in names(h5)){
      if(i == 'obs'){
        # --obs
        obs.df <- h5_to_df(h5df = h5[[i]])
        obs.name <- rownames(obs.df)
        cds_list[[i]] <- obs.df
        cds_list[[paste0(i,'.name')]] <- rownames(obs.df)
      }
      if(i == 'var'){
        # -- var
        var.df <- h5_to_df(h5df = h5[[i]])
        var.name <- rownames(var.df)
        cds_list[[i]] <- var.df
        cds_list[[paste0(i,'.name')]] <- rownames(var.df)
      }
      if(i == 'rawData'){
        rdata = h5_to_matrix(h5 = h5[[i]])
        cds_list[[i]] <- rdata
      }
      if(i == 'normData'){
        ndata = h5_to_matrix(h5 = h5[[i]])
        cds_list[[i]] <- ndata
      }
      if(i == 'dimR'){
        dimR <- h5[[i]]
        dimR_list <- list()
        for(DR in names(dimR)){
          dim_recover_ <- t(dimR[[DR]][,])
          # colnames(dim_recover_) <- paste0(DR,1:ncol(dim_recover_))
          dimR_list[[DR]] <- dim_recover_
        }
        cds_list[[i]] <- dimR_list
      }
    }
    # add the data dimnames
    exprs <- cds_list[['rawData']]
    #gr_list <- list(rawData = 'counts', normData = 'logcounts')
    rownames(exprs) <- cds_list[['var.name']]
    colnames(exprs) <- cds_list[['obs.name']]
    # -- create the monocle object: add the counts or norm data
    cds <- new_cell_data_set(expression_data = exprs,
                             cell_metadata = cds_list[['obs']],
                             gene_metadata = cds_list[['var']])
    # -- add the reduction dimension
    if('dimR' %in% names(cds_list)){
      dimR_list = cds_list[['dimR']]
      for(DR in names(dimR_list)){
        rownames(cds_list[['dimR']][[DR]]) <- cds_list[['obs.name']]
        cds@int_colData$reducedDims[[DR]] <- cds_list[['dimR']][[DR]]
      }
    }
  }
  return(cds)
}

# --- monocle write h5 file

#' The monocle is converted to h5 file
#'
#' monocle object is converted to h5 file.
#' @param sce The monocle object.
#' @param file The h5 file.
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#'
cds_write_h5 <- function(cds = NULL, file = NULL, assay.name = NULL){
  if(is.null(file)){
    stop('No such file or directory')
  }
  if(class(cds) != 'cell_data_set'){
    stop('object ', substitute(cds), ' class is not SingleCellExperiment object')
  }
  h5 <- H5File$new(filename = file, mode = 'w')
  tryCatch({
    cds_to_h5(cds = cds, h5 = h5, assay.name = assay.name)
    h5attr(h5, 'assay_name') <- assay.name
  },
  error = function(e){
    print(e)
    file.remove(file)
  },
  finally = {
    h5$colse_all()
  })
}

#' The monocle is converted to h5 file
#'
#' monocle object is converted to h5 file.
#' @param sce The monocle object.
#' @param h5 The h5 file in R
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#'
cds_to_h5 <- function(cds, h5, assay.names='RNA'){
  h5attr(h5, 'assay_name') <- assay.name
  #--- save matrix
  matrix_to_h5(mat = exprs(cds), h5 = h5, gr_name = 'rawData',save.obs.name = FALSE, save.var.name = FALSE)
  #--- save cell annotation
  df_to_h5(df = pData(cds), h5 = h5, gr_name = 'obs')
  #--- save gene annotattion
  df_to_h5(df = fData(cds), h5 = h5, gr_name = 'var')
  #--- save reduction dimension
  if(length(cds@int_colData$reducedDims)>0){
    dimReduction <- h5$create_group('dimR')
    for(D in names(cds@int_colData$reducedDims)){
      dimReduction[[D]] <- cds@int_colData$reducedDims[[D]]
    }
  }
  # to be continues
}

