#!/usr/bin/env Rscript


option_list = list(

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="RAW rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--output"), type="character", default="output.rds", 
              help="Output RDS file name", metavar="character")

)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)


Spatial_Data=readRDS(opt$rds)


SCTransform(Spatial_Data,assay = "Spatial") -> Spatial_Data


DefaultAssay(Spatial_Data) <- "Spatial"

saveRDS(Spatial_Data,file = opt$output)