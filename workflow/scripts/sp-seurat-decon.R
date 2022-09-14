#!/usr/bin/env Rscript
option_list = list(


    optparse::make_option(c("--sprds"), type="character", default=NULL, 
              help="Processed rds file of a spatial object", metavar="character"),

    optparse::make_option(c("--scrds"), type="character", default=NULL, 
              help="Processed rds file of a single cell object", metavar="character"),
        optparse::make_option(c("--output"), type="character", default="seuratdecon.pdf", 
              help="Output pdf file name", metavar="character"),
                  optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character")



)


source("workflow/scripts/scrna-functions.R")

opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$sprds) || is.null(opt$scrds) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(Seurat)
require(tidyverse)
require(viridis)

params=list(k.anchor=23,k.score=5,k.filter=100,n.trees=100)



Spatial_Data=readRDS(opt$sprds)
scrna_data=readRDS(opt$scrds)


tryCatch(scrna_data[["SCT"]],error=function(e) {

SCTransform(scrna_data) -> scrna_data 

return(scrna_data)

}) -> scrna_data

DefaultAssay(scrna_data) <- "SCT"


function_image_fixer(Spatial_Data,opt$sampleid) -> Spatial_Data


function_decon_seurat = function(reference,query,anc){

anchors <- FindTransferAnchors(reference = reference, query = query, normalization.method = "SCT",
k.anchor = anc,
k.score = params$k.score,
k.filter=params$k.filter,
n.trees = params$n.trees,recompute.residuals = FALSE)

predictions.assay <- TransferData(anchorset = anchors, refdata = reference$seurat_clusters, prediction.assay = TRUE)
query[["predictions"]] <- predictions.assay


return(query)
}



for (i in c(params$k.anchor,20,15,10,5)) {

try({
  function_decon_seurat(reference=scrna_data,query=Spatial_Data,anc=i) -> Spatial_Data
  break
  }
  )


}

DefaultAssay(Spatial_Data) <- "predictions"

cell_types_all=Idents(scrna_data) %>% unique() %>% as.character()

wp=seurat_plotting()

ggsave(opt$output,wp,height=18,width=8)
