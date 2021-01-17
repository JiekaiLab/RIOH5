# --- SingleCellExperiment(sce) H5 IO

# library
library(hdf5r)
library(Matrix)
library(SingleCellExperiment)

# --- SingleCellExperiment write h5 file

#' The singlecellexperiment is converted to h5 file
#'
#' Singlecellexperiment object is converted to h5 file.
#' @param sce The singlecellexperiment object.
#' @param file The h5 file.
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#'
sce_write_h5 <- function(sce = NULL, file = NULL, assay.name = NULL){
  if(is.null(file)){
    stop('No such file or directory')
  }
  if(class(sce) != 'SingleCellExperiment'){
    stop('object ', substitute(sce), ' class is not SingleCellExperiment object')
  }
  h5 <- H5File$new(filename = file, mode = 'w')
  tryCatch({
    sce_to_h5(sce = sce, h5 = h5, assay.name = assay.name)
    h5attr(h5, 'assay_name') <- assay.name
  },
  error = function(e){
    print(e)
  },
  finally = {
    h5$close_all()
  })
}

#' The singlecellexperiment is converted to h5 file
#'
#' Singlecellexperiment object is converted to h5 file.
#' @param sce The singlecellexperiment object.
#' @param h5 The h5 file in R.
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#'
sce_to_h5 <- function(sce, h5, assay.name){
  h5attr(h5, 'assay_name') <- assay.name
  #--- save the matrix
  gr_list <- list(counts = 'rawData', logcounts = 'normData')
  for(i in names(gr_list)){
    matrix_to_h5(mat = assay(sce, i), h5 = h5, gr_name = gr_list[[i]])
  }
  #--- save cell annotation
  df_to_h5(df = colData(sce), h5 = h5, gr_name = 'obs')
  #--- save var
  df_to_h5(df = rowData(sce), h5 = h5, gr_name = 'var')
  #--- save reduction dimension
  if(length(reducedDimNames(sce))>0){
    dimReduction <- h5$create_group('dimR')
    for(d in reducedDimNames(sce)){
      D = toupper(d)
      dimReduction[[D]] <- t(reducedDim(sce, d))
    }
  }

}

# --- singlecellexperiment read the h5 file

#' H5 to singlecellexperiment object
#'
#' Read h5 and converted h5 to the singlecellexperiment object
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#'
sce_read_h5 <- function(file=NULL, assay.name = NULL){
  if(!file.exists(file)){
    stop('No such file or directory')
  }
  h5 <- H5File$new(filename = file, mode = 'r')
  tryCatch({
    data <- h5_to_sce(h5 = h5, assay.name = assay.name)
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


#' H5 to the singlecellexperiment
#'
#' H5 file is converted to the singlecellexperiment obejct
#' @param h5 The h5 file in R
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
h5_to_sce <- function(h5, assay.name ='RNA'){
  options (warn = -1)
  if(h5attr(h5, 'assay_name') == assay.name){
    #--- create sce module list
    sce_list <- list()
    for(i in names(h5)){
      if(i == 'obs'){
        #--- read the cell annotation
        obs.df <- h5_to_df(h5df = h5[[i]])
        obs.name <- rownames(obs.df)
        sce_list[[i]] <- obs.df
        sce_list[[paste0(i,'.name')]] <- rownames(obs.df)
      }
      if(i == 'var'){
        #--- read the gene annotation
        var.df <- h5_to_df(h5df = h5[[i]])
        var.name <- rownames(var.df)
        sce_list[[i]] <- var.df
        sce_list[[paste0(i,'.name')]] <- rownames(var.df)
      }
      if(i == 'rawData'){
        #--- read the raw counts
        rdata = h5_to_matrix(h5mat = h5[[i]])
        sce_list[[i]] <- rdata
      }
      if(i == 'normData'){
        #--- read the norm data
        ndata = h5_to_matrix(h5mat = h5[[i]])
        sce_list[[i]] <- ndata
      }
      if(i == 'dimR'){
        #--- read the dimension reduction
        dimR <- h5[[i]]
        dimR_list <- list()
        for(DR in names(dimR)){
          dr <- tolower(DR)
          dim_recover_ <- t(dimR[[DR]][,])
          colnames(dim_recover_) <- paste0(DR, '_',1:ncol(dim_recover_))
          dimR_list[[dr]] <- dim_recover_
        }
        sce_list[[i]] <- dimR_list
      }
      ### to be continues
    }
    if('dimR' %in% names(sce_list)){
      #--- add the dimR names
      dimR_list = sce_list[['dimR']]
      for(dr in names(dimR_list)){
        rownames(sce_list[['dimR']][[dr]]) <- sce_list[['obs.name']]
      }
    }
    #--- add the data names
    assay_list =list()
    gr_list <- list(rawData = 'counts', normData = 'logcounts')
    for(d in grep('Data', names(sce_list), value = TRUE)){
      rownames(sce_list[[d]]) <- sce_list[['var.name']]
      colnames(sce_list[[d]]) <- sce_list[['obs.name']]
      assay_list[[gr_list[[d]]]] <- sce_list[[d]]
    }
    #--- create object
    sce = SingleCellExperiment::SingleCellExperiment(assays = assay_list,
                                                     colData = sce_list[['obs']],
                                                     rowData = sce_list[['var']],
                                                     reducedDims = sce_list[['dimR']])
  }
  return(sce)
}








