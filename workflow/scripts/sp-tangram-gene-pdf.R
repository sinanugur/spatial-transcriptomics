#!/usr/bin/env Rscript
option_list = list(

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a spatial Seurat object", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--predicted"), type="character", default=NULL, 
              help="CSV table of predictions", metavar="character"),
    optparse::make_option(c("--measured"), type="character", default=NULL, 
              help="CSV table of measured", metavar="character"),
    optparse::make_option(c("--output.dir"), type="character", default="./", 
              help="Output directory", metavar="character")


)
 
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) || is.null(opt$predicted) || is.null(opt$measured)  ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(viridis)
require(patchwork)

source("workflow/scripts/scrna-functions.R")


dir.create(opt$output.dir)
Spatial_Data=readRDS(opt$rds)
tangram_csv_predicted=read.csv(opt$predicted)
tangram_csv_measured=read.csv(opt$measured)

function_image_fixer(Spatial_Data,opt$sampleid) -> Spatial_Data



Spatial_Data[["predictions"]] <- CreateAssayObject(tangram_csv_predicted  %>% column_to_rownames("X") %>% t())
Spatial_Data[["measured"]] <- CreateAssayObject(tangram_csv_measured  %>% column_to_rownames("X") %>% t())


#DefaultAssay(Spatial_Data) <- "predictions"


cell_types_all=tangram_csv_predicted  %>% column_to_rownames("X") %>% colnames()


options(warn=-1)
suppressMessages(for (i in cell_types_all) {


try({

DefaultAssay(Spatial_Data) <- "predictions"
p1 <- SpatialFeaturePlot(Spatial_Data,features = i, ncol = 1, alpha = c(0.1, 1),images=paste0("image"),pt.size.factor=1.1) + scale_fill_continuous(type="viridis")

DefaultAssay(Spatial_Data) <- "measured"
p2 <- SpatialFeaturePlot(Spatial_Data, features = i, ncol = 1, alpha = c(0.1, 1),images=paste0("image"),pt.size.factor=1.1) + scale_fill_continuous(type="viridis")

suppressWarnings(((p2|p1)) -> wp)


ggsave(paste0(opt$output.dir,"/",i,".pdf"),wp,height=4,width=8)

})

})
