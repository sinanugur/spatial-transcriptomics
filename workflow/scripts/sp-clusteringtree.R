#!/usr/bin/env Rscript


option_list = list(
  optparse::make_option(c("--scale.factor"), type="integer", default=10000, 
              help="Scale factor [default= %default]", metavar="character"),
  optparse::make_option(c("--nfeatures"), type="integer", default=2000, 
              help="Highly variable features [default= %default]", metavar="integer"),
    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="RAW rds file of a Seurat object", metavar="character"),

    optparse::make_option(c("--normalization.method"), type="character", default="SCT", 
              help="Normalization method[default= %default]", metavar="character"),
          optparse::make_option(c("--output"), type="character", default="clustree.pdf", 
              help="Output clustree file name", metavar="character"),
                            optparse::make_option(c("--hvfplot"), type="character", default="variable-features.pdf", 
              help="Variable features file name", metavar="character"),
                                optparse::make_option(c("--heatmap"), type="character", default="dimheatmap.pdf", 
              help="Dim heatmap plot file name", metavar="character")


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(tidyverse)
require(optparse)
require(Seurat)
require(clustree)


source("workflow/scripts/scrna-functions.R")

scrna=readRDS(file = opt$rds)

scrna <- SCTransform(scrna,assay = "Spatial",verbose = FALSE)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = opt$nfeatures)


scrna <- RunPCA(scrna,assay="SCT",features = VariableFeatures(object = scrna))



#output.dir=paste0("results/",opt$sampleid,"/technicals/")
#dir.create(output.dir,recursive = T)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scrna), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scrna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)



#ggsave(paste0(output.dir,"highly-variable-features.pdf"), plot2 ,width = 8,height = 9)
ggsave(opt$hvfplot, plot2 ,width = 8,height = 9)

DimHeatmap(scrna, dims = 1:15, cells = 500, balanced = TRUE,fast = FALSE)
#ggsave(paste0(output.dir,"DimHeatMap_plot.pdf") ,width = 8,height = 15)
ggsave(opt$heatmap ,width = 8,height = 15)


dimensionReduction=function_pca_dimensions(scrna)
scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)
scrna <- FindClusters(scrna, resolution = seq(0.1,2.5,0.1))

clustree(scrna) -> p1

#output.dir=paste0("results/",opt$sampleid,"/clusteringTree/")
#dir.create(output.dir,recursive = T)

#ggsave(paste0(output.dir,"/clusteringTree-",opt$sampleid,".pdf"),p1,width=8,height=15)
ggsave(opt$output,p1,width=8,height=15)

