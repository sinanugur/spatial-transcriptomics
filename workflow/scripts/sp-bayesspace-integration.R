#!/usr/bin/env Rscript

option_list = list(
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("-r","--rds"), type="character", default=NULL, 
              help="A list of RDS files of Seurat objects", metavar="character"),
    optparse::make_option(c("--output"), type="character", default=NULL, 
              help="Output RDS file name", metavar="character"),
    optparse::make_option(c("--umap.plot"), type="character", default="umap.plot.pdf", 
              help="UMAP plot file name", metavar="character"),
    optparse::make_option(c("--harmony.plot"), type="character", default="harmony.plot.pdf", 
              help="Harmony pplot file name", metavar="character")


    


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds and sampleid)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(BayesSpace)
require(scater)
require(harmony)



files= unlist(strsplit(opt$rds, " "))
print(files)
for(i in files) {

    if(!exists("sce.combined")) {

        sce_ = readRDS(file = i)
        rowData(sce_)$is.HVG=NULL

        sce.combined = sce_

    } else {
        sce_ = readRDS(file = i)
        rowData(sce_)$is.HVG=NULL
        sce.combined= cbind(sce.combined,sce_,deparse.level = 1)

        
    }
    


}

#Combine into 1 SCE and preprocess
sce.combined = spatialPreprocess(sce.combined, n.PCs = 25,platform="Visium") #lognormalize, PCA
sce.combined =  runUMAP(sce.combined, dimred = "PCA")

colnames(reducedDim(sce.combined, "UMAP")) = c("UMAP1", "UMAP2")


ggplot(data.frame(reducedDim(sce.combined, "UMAP")) %>% dplyr::mutate(sample=factor(sce.combined$sample_name)), 
       aes(x = UMAP1, y = UMAP2, color = sample)) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw() + facet_wrap(~sample) -> p1

ggsave(filename = opt$umap.plot,p1)

sce.combined = RunHarmony(sce.combined, "sample_name", verbose = F)
sce.combined = runUMAP(sce.combined, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce.combined, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce.combined, "UMAP.HARMONY")) %>% dplyr::mutate(sample=factor(sce.combined$sample_name)), 
       aes(x = UMAP1, y = UMAP2, color = sample)) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw() + facet_wrap(~sample,ncol = 4) -> p2

ggsave(filename = opt$harmony.plot,p1)

saveRDS(sce.combined,file= opt$output)