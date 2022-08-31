#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--nfeatures"), type="integer", default=2000, 
              help="Highly variable features [default= %default]", metavar="integer"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="RAW rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character")


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)

source("workflow/scripts/scrna-functions.R")

scrna=readRDS(file = opt$rds)


scrna <- SCTransform(scrna,assay = "Spatial",verbose = FALSE)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = opt$nfeatures)


scrna <- RunPCA(scrna,assay="SCT", features = VariableFeatures(object = scrna))
dimensionReduction=function_pca_dimensions(scrna)
scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)
scrna <- FindClusters(scrna, resolution = opt$resolution)
scrna <- RunUMAP(scrna, dims = 1:dimensionReduction)




output.dir=paste0("analyses/processed/",opt$resolution,"/")
dir.create(output.dir,recursive = T)

saveRDS(scrna,file = paste0(output.dir,opt$sampleid,".rds"))


SCT_=paste0("SCT_snn_res.",opt$resolution)

metrics=table(scrna@meta.data[[SCT_]], scrna@meta.data$orig.ident)

output.dir=paste0("results/",opt$sampleid,"/resolution-",opt$resolution,"/")
dir.create(output.dir,recursive = T)

openxlsx::write.xlsx(metrics %>% as.data.frame() %>% select(Cluster=1,everything()),file=paste0(output.dir,opt$sampleid,".number-of-cells-per-cluster",".xlsx"))



