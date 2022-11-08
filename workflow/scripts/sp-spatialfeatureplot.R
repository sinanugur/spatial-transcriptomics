#!/usr/bin/env Rscript
option_list = list(

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--output"), type="character", default="output.pdf", 
              help="Output pdf file name", metavar="character")


)
 
 source("workflow/scripts/scrna-functions.R")



opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(viridis)
require(patchwork)


Spatial_Data=readRDS(opt$rds)


#plot1 <- VlnPlot(Spatial_Data, features = c("nCount_Spatial","nFeature_Spatial"), ncol = 2) & theme(legend.position = "bottom")
plot2 <- SpatialFeaturePlot(Spatial_Data, features = c("nCount_Spatial","nFeature_Spatial"),pt.size.factor=4.3,ncol = 2) & 
theme(legend.position = "right") & scale_fill_continuous(type="viridis")
#wrap_plots(plot1, plot2) -> wp


ggsave(opt$output,plot2,height=4,width=8.5)