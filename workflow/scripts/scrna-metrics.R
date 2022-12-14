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



source("workflow/scripts/scrna-functions.R")


scrna=readRDS(file = opt$rds)


output.dir=paste0("results/",opt$sampleid,"/technicals/")
dir.create(output.dir,recursive = T)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scrna), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scrna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)



ggsave(paste0(output.dir,"highly-variable-features.pdf"), plot2 ,width = 8,height = 9)




DimHeatmap(scrna, dims = 1:15, cells = 500, balanced = TRUE,fast = FALSE)
ggsave(paste0(output.dir,"DimHeatMap_plot.pdf") ,width = 8,height = 15)



scrna <- JackStraw(scrna, num.replicate = 100,  dims=50)
scrna <- ScoreJackStraw(scrna, dims = 1:50)


plot1 <- JackStrawPlot(scrna, dims = 1:50) 
plot2 <- ElbowPlot(scrna, ndims=50)
ggsave(paste0(output.dir,"JackandElbow_plot.pdf"), plot1 + plot2,width = 13,height = 5)