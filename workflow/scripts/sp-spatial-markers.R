#!/usr/bin/env Rscript
option_list = list(

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="Processed rds file of a Seurat object", metavar="character"),

    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--selection.method"), type="character", default="markvariogram", 
              help="Spatial marker selection method", metavar="character"),
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


Spatial_Data=readRDS(opt$rds)

function_image_fixer(Spatial_Data,opt$sampleid) -> Spatial_Data

Spatial_Data <- SCTransform(Spatial_Data, assay = "Spatial", verbose = FALSE)
Spatial_Data <- FindSpatiallyVariableFeatures(Spatial_Data, assay = "SCT", features = VariableFeatures(Spatial_Data)[1:1000],
    selection.method = opt$selection.method)

markers=SpatiallyVariableFeatures(Spatial_Data, selection.method = opt$selection.method) 

openxlsx::write.xlsx(markers %>% as.data.frame() %>% select(gene=1) %>% head(50),file=opt$output)

output.dir=paste0("results/",opt$sampleid,"/spatial-markers/plots/")
dir.create(output.dir,recursive = T)

options(warn=-1)
for (i in markers) {

#SpatialFeaturePlot(Spatial_Data, features = i, ncol = 1, alpha = c(0.1, 1),images=paste0("image"),pt.size.factor=1.1) + scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn")))

SpatialFeaturePlot(Spatial_Data, features = i, ncol = 1, alpha = c(0.1, 1),images=paste0("image"),pt.size.factor=1.1) + scale_fill_continuous(type="viridis")
ggsave(paste0("results/",opt$sampleid,"/","spatial-markers/plots/",i,".pdf"))

}