#!/usr/bin/env Rscript
option_list = list(

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),

    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
      optparse::make_option(c("--n"), type="integer", default=50, 
              help="How many features to plot per cluster [default= %default]", metavar="integer"),
    optparse::make_option(c("--csv"), type="character", default=NULL, 
              help="Excel table of markers", metavar="character"),
    optparse::make_option(c("--output"), type="character", default="tangram.pdf", 
              help="Output tangram file name", metavar="character")


)
 
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) || is.null(opt$csv)  ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(viridis)


source("workflow/scripts/scrna-functions.R")



Spatial_Data=readRDS(opt$rds)
tangram_csv=read.csv(opt$csv)



function_image_fixer(Spatial_Data,opt$sampleid) -> Spatial_Data

Spatial_Data[["predictions"]] <- CreateAssayObject(tangram_csv  %>% column_to_rownames("X") %>% t())


DefaultAssay(Spatial_Data) <- "predictions"


cell_types_all=tangram_csv  %>% column_to_rownames("X") %>% colnames()


wp=seurat_plotting()


ggsave(opt$output,wp,height=18,width=8)
