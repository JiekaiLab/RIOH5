# Seurat object H5 IO

# library
library(Seurat)
library(hdf5r)
library(Matrix)

# --- seurat read the H5 file

#' H5 to Seuart object
#'
#' Read h5 and converted h5 to the seurat object
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#'
seurat_read_h5 <- function(file=NULL, assay.name = NULL){
  if(!file.exists(file)){
    stop('No such file or directory')
  }
  h5 <- H5File$new(filename = file, mode = 'r')
  tryCatch({
    data <- h5_to_seurat(h5 = h5, assay.name = assay.name)
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


#' H5 to the seurat
#'
#' H5 file is converted to the seurat obejct
#' @param h5 The h5 file in R
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @param load.graph The knn and snn network isn't loaded(FALSE by defualt).
#'
h5_to_seurat <- function(h5, assay.name, load.graph=FALSE){
  options(warn = -1)
  if(h5attr(h5, 'assay_name') == assay.name){
    # -- create seurat module list
    seurat_list <- list()
    for(i in names(h5)){
      if(i == 'obs'){
        #--- read cell annotation
        # obs.df <- h5_to_df(h5 = h5[[i]])
        # obs.name <- rownames(obs.df)
        seurat_list[[i]] <- h5_to_df(h5 = h5[[i]])
        seurat_list[[paste0(i,'.name')]] <- h5[[i]][['index']][]
      }
      if(i == 'var'){
        #--- read gene annotation
        # var.df <- h5_to_df(h5 = h5[[i]])
        # var.name <- rownames(var.df)
        seurat_list[[i]] <- h5_to_df(h5 = h5[[i]])
        seurat_list[[paste0(i,'.name')]] <- h5[[i]][['index']][]
      }
      if(i == 'rawData'){
        #--- read the raw counts
        # rdata = h5_to_matrix(h5 = h5[[i]])
        seurat_list[[i]] <- h5_to_matrix(h5 = h5[[i]])
      }
      if(i == 'normData'){
        #--- read the norm data
        # ndata = h5_to_matrix(h5 = h5[[i]])
        seurat_list[[i]] <- h5_to_matrix(h5 = h5[[i]])
      }
      if(i == 'dimR'){
        #--- read the dimension reduction
        dimR <- h5[[i]]
        dimR_list <- list()
        for(DR in names(dimR)){
          dr <- tolower(DR)
          # dim_recover <- t(dimR[[DR]][,])
          dim_recover_ <- Seurat::CreateDimReducObject(embeddings = t(dimR[[DR]][,]), key = paste0(DR, "_"), assay = "RNA")
          dimR_list[[dr]] <- dim_recover_
        }
        seurat_list[[i]] <- dimR_list
      }
      ### to be continues
    }
  }else{stop("Entering the 'assay.name' corresponding to h5")}
  #--- add the data names
  for(d in grep('Data', names(seurat_list), value = TRUE)){
    rownames(seurat_list[[d]]) <- seurat_list[['var.name']]
    colnames(seurat_list[[d]]) <- seurat_list[['obs.name']]
  }
  #--- create the seurat object: add the counts or norm data
  if(all(c('normData', 'rawData')%in%names(seurat_list))){
    seurat <- Seurat::CreateSeuratObject(counts = seurat_list[['rawData']], assay = assay.name)
    seurat@assays[[assay.name]]@data <- seurat_list[['normData']]
  }
  else if(sum(c('normData', 'rawData')%in%names(seurat_list)) == 1){
    if('rawData' %in% names(seurat_list)){
      seurat <- Seurat::CreateSeuratObject(counts = seurat_list[['rawData']], assay = assay.name)
    }
    else if('normData'%in% names(seurat_list)){
      seurat <- Seurat::CreateSeuratObject(counts = seurat_list[['normData']], assay = assay.name)
    }
  }
  #--- create the seurat object: add the meta.data and meta.features
  slot(object = seurat, name = 'meta.data') <- seurat_list[['obs']]
  seurat@assays[[assay.name]]@meta.features <- seurat_list[['var']]
  #--- create the seurat object: add the dimension reduction
  if('dimR' %in% names(seurat_list)){
    for(dr in names(seurat_list[['dimR']])){
      seurat@reductions[[dr]] <- seurat_list[['dimR']][[dr]]
      rownames(seurat@reductions[[dr]]@cell.embeddings) = rownames(seurat[[]])
    }
  }
  return(seurat)
}

# --- seurat write the h5 file

#' The seurat object to h5
#'
#' The Seurat object is converted to the h5 file.
#' @param seurat The seurat object.
#' @param file The h5 file.
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @param save.graphs The knn and snn network isn't saved(FALSE by defualt).
#' @param save.scale.data The scale data isn't saved(FALSE by defualt).
#'
seurat_write_h5 <- function(seurat = NULL, file = NULL, assay.name = NULL, save.graphs = FALSE, save.scale.data = FALSE){
  if(is.null(file)){
    stop('No such file or directory')
  }
  if(class(seurat) != 'Seurat'){
    stop('object ', substitute(seurat), ' class is not Seurat object')
  }
  h5 <- H5File$new(filename = file, mode = 'w')
  tryCatch({
    seurat_to_h5(seurat = seurat, h5 = h5, assay.name = assay.name, save.graphs = save.graphs, save.scale.data = save.scale.data)
    h5attr(h5, 'assay_name') <- assay.name
  },
  error = function(e) {
    print(e)
  },
  finally = {
    h5$close_all()
  }
  )
}


#' Seurat is converted to h5 file
#'
#' Seurat object is converted to h5 file.
#' @param seurat The seurat object.
#' @param h5 The h5 file in R.
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @param save.graphs The knn and snn network isn't saved(FALSE by defualt).
#' @param save.scale.data The scale data isn't saved(FALSE by defualt).
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#'
seurat_to_h5 <- function(seurat=NULL, h5=NULL, assay.name = NULL, save.graphs = FALSE, save.scale.data = FALSE){
  if(assay.name %in% names(slot(object = seurat, name = 'assays'))){
    #--- data
    slot_assay <- slot(object = seurat, name = 'assays')[[assay.name]]
    if(all(slot_assay@counts@x == slot_assay@data@x)){
      # only have rdata (rdata == ndata) defualt norm data
      rdata <- slot(object = slot_assay, name = 'counts')
      matrix_to_h5(mat = rdata, h5 = h5, gr_name = 'rawData')
      print("'counts' is the sames as the 'data'.")
    }
    else{
      # have rdata and ndata (rdata != ndata)
      rdata <- slot(object = slot_assay, name = 'counts')
      matrix_to_h5(mat = rdata, h5 = h5, gr_name= 'rawData')
      ndata <- slot(object = slot_assay, name = 'data')
      matrix_to_h5(mat = ndata, h5 = h5, gr_name = 'normData')
    }
    #--- save the cell annotation
    obs = slot(object = seurat, name = 'meta.data')
    df_to_h5(df = obs, h5 = h5, gr_name = 'obs')
    #--- save the gene annotation
    ndvar <- slot(object = slot_assay, name = 'meta.features')
    hvg <- slot(object = slot_assay, name = 'var.features')
    ndvar[['highly_variable']] <- rownames(ndvar) %in% hvg
    df_to_h5(df = ndvar, h5 = h5, gr_name= 'var')
    #--- save the dimension reduction
    if(length(seurat@reductions)>0){
      dimReduction <- h5$create_group('dimR')
      dimr <- slot(object = seurat, name = 'reductions')
      for(d in names(dimr)){
        D = toupper(d)
        dimReduction[[D]] <- t(slot(dimr[[d]],'cell.embeddings'))
      }
    }
    #--- save the graphs
    if(save.graphs){
      if(length(seurat@graphs)>0){
        graph_df <- slot(object = seurat, 'graphs')
        if(length(grep(assay.name, names(graph_df))) == 2){
          graphs <- h5$create_group('graphs')
          gra_list <- list()
          gra_list[[paste0(assay.name, '_nn')]] <- 'knn'
          gra_list[[paste0(assay.name, '_snn')]] <- 'snn'
          for(g in names(graph_df)){
            matrix_to_h5(mat=graph_df[[g]], h5 = graphs, gr_name = gra_list[[g]])
          }
        }
      }
    }
  }
  else{
    stop('Please enter the correct assay name')
  }
  return('Done!')
}






