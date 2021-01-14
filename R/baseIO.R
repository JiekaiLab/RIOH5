# ---  base data IO

# packages
library(hdf5r)
library(Matrix)

#' Data frame to h5
#'
#' Data frame is converted to and saved the h5 file
#' @param df Data frame of cell annotation or gene annotation
#' @param h5 The h5 file
#' @param gr_name The group name represents the property of the data frame.
#' @export
#'
df_to_h5 <- function(df, h5, gr_name=NULL){
  # if(!gr_name %in% names(h5)){
  #   h5df <- h5$create_group(gr_name)
  # }
  # else{h5df <- h5[[gr_name]]
  # }
  h5df <- h5$create_group(gr_name)
  h5df[['index']] = rownames(df)
  h5df[['colnames']] = colnames(df)
  # factor to levels,charactor to levels,logical to levels
  for(k in names(df)){
    if(is.factor(df[[k]])){
      h5df[[k]]<- as.integer(df[[k]]) - 1L # for 0 begin
      h5df[[paste0(k,'_levels')]]<- levels(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'category'
    }
    if(is.character(df[[k]])){
      str_to_lvl <- factor(df[[k]])
      h5df[[k]]<- as.integer(str_to_lvl) - 1L
      h5df[[paste0(k,'_levels')]]<- levels(str_to_lvl)
      h5attr(h5df[[k]], 'origin_dtype') = 'string'
    }
    if(is.logical(df[[k]])){
      h5df[[k]] <- as.integer(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'bool'
    }
    if(any(is.numeric(df[[k]]),is.integer(df[[k]]))){
      h5df[[k]] <- df[[k]]
      h5attr(h5df[[k]], 'origin_dtype') = 'number'
    }
  }
}

#' H5 to dataframe
#'
#' H5 gruop is converted to usable data frame including the observes annotation (cell annotation) and variables annotation (gene annotation)
#' @param h5df The hdf5 group saved the data frame using list mode.
#' @return A data frame(annotation)
#' @export
#'
h5_to_df <- function(h5df){
  df_list <- list()
  df_list[['index']] <- h5df[['index']][]
  cnames <- h5df[['colnames']][]
  for(k in names(h5df)){
    if(length(h5attr_names(h5df[[k]]))>0){
      df_dtype <- h5attr(h5df[[k]], 'origin_dtype')
      if(df_dtype == 'category'){
        e0 <- h5df[[k]][] + 1L
        lvl <- h5df[[paste0(k,'_levels')]][]
        df_list[[k]] <- structure(.Data = e0, .Label = lvl, class = 'factor')
      }
      if(df_dtype == 'string'){
        e0 <- h5df[[k]][] + 1L
        lvl <- h5df[[paste0(k,'_levels')]][]
        df_list[[k]] <- as.character(structure(.Data = e0, .Label = lvl, class = 'factor'))
      }
      if(df_dtype == 'bool'){
        df_list[[k]] <- as.logical(h5df[[k]][])
      }
      if(df_dtype == 'number'){
        df_list[[k]] <- h5df[[k]][]
      }
    }
  }
  df = as.data.frame(df_list, row.names = 'index')
  if(length(cnames)>1){
    df = df[,cnames]
  }else{
    df=df
  }
  return(df)
}


#' Matrix to H5 format
#'
#' The matrix including the dense matrix and sparse matrix is converted to the matrix in h5 format or is stored into the h5 file.
#' @param mat The matrix object including matrix(R) and sparse matrix(Matrix package)
#' @param h5 The H5 file name that we write in
#' @param gr_name The h5 gorup store the matrix (dense matrix or sparse matrix)
#' @param save.obs.name The rownames isn't be saved(FALSE by defualt)
#' @param save.var.name The colnames isn't be saved(FALSE by defualt)
#' @export
#'
matrix_to_h5 <- function(mat, h5, gr_name = NULL, save.obs.name = FALSE, save.var.name = FALSE){
  if(!gr_name %in% names(h5)){
    h5mat = h5$create_group(gr_name)
  }
  else{
    h5mat = h5[[gr_name]]
  }
  if('dgCMatrix' %in% class(mat)){
    h5mat[['values']] <- slot(object = mat, name = 'x')
    h5mat[['indices']] <- slot(object = mat, name = 'i')
    h5mat[['indptr']] <- slot(object = mat, name = 'p')
    h5mat[['dims']] <- rev(slot(object = mat, name = 'Dim'))
    if(save.obs.name & save.var.name){
      h5mat[['var_names']] <- slot(object = mat, name = 'Dimnames')[[1]]
      h5mat[['obs_names']] <- slot(object = mat, name = 'Dimnames')[[2]]
    }
    h5attr(h5mat, 'datatype') <- 'SparseMatrix'
    print(paste0(substitute(gr_name),' is sparse matrix'))#1
  }
  else if('matrix' %in% class(mat)){
    h5mat[['matrix']] <- mat
    h5mat[['dims']] <- rev(dim(mat))
    if(save.obs.name & save.var.name){
      h5mat[['var_names']] <- slot(object = mat, name = 'Dimnames')[[1]]
      h5mat[['obs_names']] <- slot(object = mat, name = 'Dimnames')[[2]]
    }
    h5attr(h5mat, 'datatype') <- 'Array'
    warning(paste0(substitute(gr_name), ' is dense matrix'))#2
  }
  else if('Graph' %in% class(mat)){
    h5mat[['values']] <- slot(object = mat, name = 'x')
    h5mat[['indices']] <- slot(object = mat, name = 'i')
    h5mat[['indptr']] <- slot(object = mat, name = 'p')
    h5mat[['dims']] <- rev(slot(object = mat, name = 'Dim'))
    h5attr(h5mat, 'datatype') <- 'SparseMatrix'
  }
}

#' H5 to Matrix format
#'
#' H5 group is converted to usable matrice including dense matrices and sparse matrices.
#' Dense matrices(by R) and sparse matrices(constructed by Matrix package)
#' @param h5 The name of the group the stores the matrix data in h5 file
#' @param obs_names The observes names, such as cell names
#' @param var_names The variables names, such as gene names
#' @return dense matrix or sparse matrix
#' @export
#'
h5_to_matrix <-  function(h5mat, obs.name=NULL, var.name=NULL){
  if(all(c('obs_names','var_names') %in% names(h5mat))){
    obs.name = h5mat[['obs_names']][]
    var.name = h5mat[['var_names']][]
  }
  if(h5attr(h5mat, 'datatype') == 'SparseMatrix'){
    mat <- Matrix::sparseMatrix(i = h5mat[['indices']][],
                                p = h5mat[['indptr']][],
                                x = h5mat[['values']][],
                                dims = rev(h5mat[['dims']][]),
                                index1 = FALSE)
    if(!is.null(obs.name) & !is.null(var.name)){
      dimnames(mat) <- list(var.name, obs.name)
    }
    else{
      warning('There are no dimnames in the sparse matrix')
    }
  }
  else if(h5attr(h5mat, 'datatype') == 'Array'){
    mat <- h5mat[['matrix']][,]
    if(!is.null(obs.name) & !is.null(var.name)){
      dimnames(mat) <- list(var.name, obs.name)
    }
    else{
      warning('There are no dimnames in the matrix')
    }
  }
  return(mat)
}

