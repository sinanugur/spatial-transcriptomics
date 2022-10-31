#!/usr/bin/env Rscript


option_list = list(

    optparse::make_option(c("--rds"), type="character", default=NULL, 
              help="RAW rds file of a Seurat object", metavar="character"),
    optparse::make_option(c("--output"), type="character", default="output.rds", 
              help="Output RDS file name", metavar="character"),
    optparse::make_option(c("--assay"), type="character", default="RNA", 
              help="Assay type, RNA or Spatial", metavar="character"),
    optparse::make_option(c("--output.xlsx"), type="character", default="output.xlsx", 
              help="Output Excel file name for cluster markers", metavar="character")



)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$rds) ){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)


if(opt$assay == "Spatial") {

Spatial_Data=readRDS(opt$rds)

CreateSeuratObject(Spatial_Data[["Spatial"]],assay="Spatial") -> Spatial_Data_tmp

Spatial_Data_tmp@images <- Spatial_Data@images

SCTransform(Spatial_Data_tmp,assay = "Spatial") -> Spatial_Data_tmp

DefaultAssay(Spatial_Data_tmp) <- "SCT"

saveRDS(Spatial_Data_tmp,file = opt$output)

} else {

scrna_data=readRDS(opt$rds)


SCTransform(scrna_data,ncells = 3000, verbose = FALSE) -> scrna_data 


DefaultAssay(scrna_data) <- "SCT"

saveRDS(scrna_data,file=opt$output)



cluster_markers_all <- Seurat::FindAllMarkers(object = scrna_data, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)

openxlsx::write.xlsx(cluster_markers_all,file=opt$output.xlsx)

}
 
