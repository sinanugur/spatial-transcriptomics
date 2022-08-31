#!/usr/bin/env Rscript
option_list = list(

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--output"), type="character", default="output.xlsx", 
              help="Output excel file name", metavar="character")


)
 
 source("workflow/scripts/scrna-functions.R")



opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(viridis)
require(patchwork)


Spatial_Data=readRDS(opt$rds)

function_image_fixer(Spatial_Data,opt$sampleid) -> Spatial_Data



plot1 <- VlnPlot(Spatial_Data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Spatial_Data, features = "nCount_Spatial",images=paste0("image"),pt.size.factor=1.1) + theme(legend.position = "right",legend.title = element_blank())
wrap_plots(plot1, plot2) -> wp


ggsave(opt$output,wp,height=4.5,width=7)