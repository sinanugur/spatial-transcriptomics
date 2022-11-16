#!/usr/bin/env Rscript

option_list = list(
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character"),
    optparse::make_option(c("-r","--rds"), type="character", default=NULL, 
              help="A list of RDS files of Seurat objects", metavar="character"),
    optparse::make_option(c("--resolution"), type="double", default=0.8, 
              help="Resolution [default= %default]", metavar="character"),
    optparse::make_option(c("--umap.plot"), type="character", default="umap.pdf", 
              help="Resolution [default= %default]", metavar="character"),
    optparse::make_option(c("--output"), type="character", default="output.pdf", 
              help="Output file name [default= %default]", metavar="character")
    


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) || is.null(opt$sampleid) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds and sampleid)", call.=FALSE)
}

require(tidyverse)
require(Seurat)
require(patchwork)
require(randomcoloR)
source("workflow/scripts/scrna-functions.R")

files= unlist(strsplit(opt$rds, " "))
print(files)
for(i in files) {

    if(!exists("scrna_list")) {

        scrna_list=list(readRDS(file = i))

    } else {

        scrna_list=append(scrna_list,readRDS(file = i))
    }
    


}


for (i in 1:length(scrna_list)){

    DefaultAssay(scrna_list[[i]]) <- "SCT"

}

features <- SelectIntegrationFeatures(object.list = scrna_list, nfeatures = 3000)
scrna_list <- PrepSCTIntegration(object.list = scrna_list, anchor.features = features)


scrna.anchors <- FindIntegrationAnchors(object.list = scrna_list, normalization.method = "SCT",
    anchor.features = features)
scrna.combined.sct <- IntegrateData(anchorset = scrna.anchors, normalization.method = "SCT")



scrna.combined.sct <- RunPCA(scrna.combined.sct, verbose = FALSE)
scrna.combined.sct  <- FindNeighbors(scrna.combined.sct , dims = 1:30)
scrna.combined.sct <- FindClusters(scrna.combined.sct , verbose = FALSE)
scrna.combined.sct <- RunUMAP(scrna.combined.sct, reduction = "pca", dims = 1:30)

palette <- distinctColorPalette(length(unique(Idents(scrna.combined.sct ))))
names(palette)=as.character(unique(Idents(scrna.combined.sct)))

p1 <- DimPlot(scrna.combined.sct, reduction = "umap", group.by = "orig.ident",cols=palette)
p2 <- DimPlot(scrna.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    repel = TRUE,cols=palette)



#output.dir=paste0("results/integration/seurat/technicals/")
#dir.create(output.dir,recursive = T)

ggsave(file=opt$umap.plot,p1+p2,height=5,width=11)




#output.dir=paste0("analyses/integration/seurat/")
#dir.create(output.dir,recursive = T)

saveRDS(scrna.combined.sct,file = opt$output)
