library(Seurat)
library(hdf5r)
library(Matrix)
library(Hmisc)


#' create the VisiumV1 class inherting the package Seurat 'SpatialImage' class
#' @import Seurat
VisiumV1 <- setClass(
  Class = 'VisiumV1',
  contains = 'SpatialImage',
  slots = list(
    'image' = 'array',
    'scale.factors' = 'scalefactors',
    'coordinates' = 'data.frame',
    'spot.radius' = 'numeric'
  )
)


#' The spatial meassage to the h5 file(Seurat)
#'
#' @param data The seruat object
#' @param h5 The h5 file
#' @param assay.name 'assay.name' must be the 'spatial' for saving the spatial message.
#'
seurat_spatial_to_h5 <- function(data, h5, gr_name = 'spatial'){
  h5spa <- h5$create_group(gr_name)
  for(sample_id in names(data@images)){
    sid_h5 <-  h5spa$create_group(sample_id)
    #--- save the image for the lowres
    sid_image_h5 <- sid_h5$create_group('image')
    sid_image_h5[['lowres']] <-  aperm(slot(data@images[[sample_id]], 'image'))
    #--- save the scale factor
    sid_scalefactors_h5 <- sid_h5$create_group('scalefactors')
    sf <- slot(data@images[[sample_id]], 'scale.factors')
    v1 <- c('spot_diameter_fullres','fiducial_diameter_fullres','tissue_hires_scalef', 'tissue_lowres_scalef')
    for(k in names(sf)){
      sid_scalefactors_h5[[grep(k,v1, value = TRUE)]] <- sf[[k]]
    }
    #--- save the coor
    coor_df = slot(data@images[[sample_id]], 'coordinates')
    colnames(coor_df) <- c('in_tissue', 'array_row','array_col', 'image_2','image_1')
    df_to_h5(df = coor_df[,c('in_tissue', 'array_row','array_col','image_1','image_2')], h5 = sid_h5, gr_name = 'coor')
  }
}


#' The h5 group spatial to the spatial message
#'
#' @param h5spa The h5 group for the spatial
#' @return The spatial message
#'
h5_to_spatial <- function(h5spa){
  spatial_list<- list()
  for(sid in names(h5spa)){
    spatial_sid_list <- list()
    sid_h5 <- h5spa[[sid]]
    for(me in names(sid_h5)){
      if('image' == me){
        spatial_sid_list[[me]] <- aperm(sid_h5[[me]][['lowres']][,,])
      }
      if('scalefactors' == me){
        sf_list <- list()
        v1 = c("spot","fiducial","hires","lowres")
        for(sf in v1){
          sf_list[[sf]] <- sid_h5[[me]][[grep(sf,names(sid_h5[[me]]), value = TRUE)]][]
        }
        sf_o <- Seurat::scalefactors(spot = sf_list$spot, fiducial = sf_list$fiducial, hires = sf_list$hiresk, lowres = sf_list$lowres)
        spatial_sid_list[[me]] <- sf_o
      }
      if('coor' == me){
        coor_df <- h5_to_df(sid_h5[['coor']])
        coor_df <- coor_df[,c('in_tissue','array_row','array_col','image_2','image_1')]
        colnames(coor_df)<- c('tissue','row','col','imagerow','imagecol')
        spatial_sid_list[[me]] <- coor_df
      }
    }
    unnormalized.radius <- spatial_sid_list$scalefactors$fiducial * spatial_sid_list$scalefactors$lowres
    spot.radius <- unnormalized.radius/max(dim(x = spatial_sid_list$image))
    spatial_list[[sid]] <- new(Class = "VisiumV1", image = spatial_sid_list$image, scale.factors = spatial_sid_list$scalefactors,
                                   coordinates = spatial_sid_list$coor,spot.radius = spot.radius)
  }
  return(spatial_list)
}

