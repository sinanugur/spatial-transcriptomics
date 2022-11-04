#!/usr/bin/env Rscript
option_list = list(


    optparse::make_option(c("--sprds"), type="character", default=NULL, 
              help="Processed rds file of a spatial object", metavar="character"),
    optparse::make_option(c("--scrds"), type="character", default=NULL, 
              help="Processed rds file of a single cell object", metavar="character"),
        optparse::make_option(c("--output"), type="character", default="output.csv", 
              help="Output CSV file name", metavar="character"),
    optparse::make_option(c("--input.xlsx"), type="character", 
              help="Input Excel file name for cluster markers", metavar="character")



)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$sprds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (sprds)", call.=FALSE)
}



require(Seurat)
require(tidyverse)
require(SPOTlight)


source("workflow/scripts/scrna-functions.R")

Spatial_Data=readRDS(opt$sprds)
scrna_data=readRDS(opt$scrds)



openxlsx::read.xlsx(opt$input.xlsx) -> cluster_markers_all

function_spotlight=function(Seurat_Object,scrna_rds) {
set.seed(123)

spotlight_ls <- spotlight_deconvolution(
  se_sc = scrna_rds,
  counts_spatial = Seurat_Object@assays$Spatial@counts,
  clust_vr = "seurat_clusters", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use #default is 100
  hvg = 3000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
  )
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(Seurat_Object)

decon_df <- decon_mtrx %>%
  data.frame() %>% select(-res_ss)
  #tibble::rownames_to_column("barcodes")

#Seurat_Object@meta.data <- Seurat_Object@meta.data %>% tibble::rownames_to_column("barcodes") %>% dplyr::left_join(decon_df, by = "barcodes") %>% tibble::column_to_rownames("barcodes")

#return(Seurat_Object)
return(decon_df)
}

function_spotlight(Spatial_Data,scrna_data) -> decon_df

write.csv(decon_df,file = opt$output)