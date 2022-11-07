#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character")


)
 


opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(patchwork)


source("workflow/scripts/scrna-functions.R")

Spatial_Data=readRDS(file = opt$rds)

function_image_fixer(Spatial_Data,opt$sampleid) -> Spatial_Data

p1 <- DimPlot(Spatial_Data, reduction = "umap", label = TRUE,label.size = 5) 
p2 <- SpatialDimPlot(Spatial_Data, label = TRUE, label.size = 5,images=paste0("image"),pt.size.factor=1.6)



output.dir=paste0("results/",opt$sampleid,"/resolution-",opt$resolution,"/")
dir.create(output.dir,recursive = T)

ggsave(plot =p1 + p2,filename=paste0(output.dir,opt$sampleid,".umap",".pdf"),width=15,height=6.5)

