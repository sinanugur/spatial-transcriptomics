#!/usr/bin/env Rscript

option_list = list(
  optparse::make_option(c("--pca.dimension"), type="integer", default=30, 
              help="PCA dimensions [default= %default]", metavar="integer"),
    optparse::make_option(c("--input"), type="character", default=NULL, 
              help="Bayesspace RDS file", metavar="character"),
    optparse::make_option(c("--output.enhanced"), type="character", default="output.enhanced.rds", 
              help="Output RDS enhanced file name", metavar="character"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("--n.cluster"), type="integer", default=NULL, 
              help="Number of clusters, if not given, autoselect", metavar="integer")
    


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$input)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (data.dir and sampleid)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)
require(BayesSpace)


domanska_muscularis=data.frame(gene=c("FCER1A","CDC1C","CLEC10A","CCL3L1","CCL3","CCL4L2","MT1X","MT1E","CTSL","RGS1","FOS","APOE","DNASE1L3","MMP9","LYZ","AREG","EREG","CCL20","S100A9","S100A8","EREG","FCN1","VCAN","LYZ","HSPA1A","HSPA6","HSPA1B","LYVE1","MARCO","COLEC12"),group=c(3,3,3,5,5,5,6,6,6,7,7,7,4,4,4,1,1,1,0,0,0,2,2,2,9,9,9,11,11,11))
domanska_mucosas=data.frame(gene=c("LYVE1","F13A1","FOLR2","SELENOP","APOE","SLC40A1","C1QB","DAB2","PDK4","SPP1","ACP5","CD9","FCER1A","CD1C","CLEC10A","HSPA6","DNAJB1","HSPA1B","S100A8","S100A9","S100A12","EREG","G0S2","FCN1","CCL20","IL1B","IL23A","CXCL10","CXCL9","GBP1"),group=c(12,12,12,11,11,11,10,10,10,9,9,9,8,8,8,5,5,5,0,0,0,3,3,3,2,2,2,6,6,6))

domanska_markers=bind_rows(domanska_mucosas,domanska_muscularis) %>% distinct(gene) %>% pull()


sce=readRDS(opt$input)


set.seed(149)
sce_enhanced <- spatialEnhance(sce, q=opt$n.cluster, platform="Visium", d=opt$pca.dimension,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=100,
                                    save.chain=TRUE)


set.seed(149)
palette <- distinctColorPalette(n)
names(palette)=as.character(unique(sce_enhanced$spatial.cluster))

clusterPlot(sce_enhanced,is.enhanced = T,palette=palette) + scale_x_reverse() + scale_y_reverse() -> wp


sce_enhanced <- enhanceFeatures(sce_enhanced, sce,
                                     feature_names=domanska_markers,model = "xgboost",
                                     nrounds=0)


sce_enhanced$sample_name <- opt$sampleid

saveRDS(sce_enhanced,opt$output.enhanced)

