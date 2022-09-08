#!/usr/bin/env Rscript

option_list = list(
    optparse::make_option(c("--type"), type="character", default="RNA", 
              help="Assay type, RNA, or Spatial", metavar="character"),
    optparse::make_option(c("-r","--rds"), type="character", default=NULL, 
              help="A list of RDS files of Seurat objects", metavar="character"),
    optparse::make_option(c("--output"), type="character", default="output.file", 
              help="Output h5ad file name", metavar="character")



)

opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call.=FALSE)
}

require(Seurat)
require(SeuratDisk)
require(tidyverse)


scrna=readRDS(file = opt$rds)

UpdateSeuratObject(scrna) -> scrna
DefaultAssay(scrna) <- opt$type

DietSeurat(scrna) -> scrna

scrna@meta.data %>% dplyr::mutate(dplyr::across(where(is.factor), as.character)) -> scrna@meta.data



SaveH5Seurat(scrna,paste0(opt$output,".h5Seurat"),overwrite = TRUE)
SeuratDisk::Convert(paste0(opt$output,".h5Seurat"), dest = "h5ad",overwrite = TRUE)