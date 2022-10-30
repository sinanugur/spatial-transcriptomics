#!/usr/bin/env Rscript
option_list = list(


    optparse::make_option(c("--sprds"), type="character", default=NULL, 
              help="Processed rds file of a spatial object", metavar="character"),
    optparse::make_option(c("--scrds"), type="character", default=NULL, 
              help="Processed rds file of a single cell object", metavar="character"),
        optparse::make_option(c("--output"), type="character", default="output.rds", 
              help="Output RDS file name", metavar="character")



)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$sprds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (sprds)", call.=FALSE)
}

source("workflow/scripts/scrna-functions.R")

Spatial_Data=readRDS(opt$spr)
dwls_data=readRDS(paste0("DWLS_assay/",scrnaID,"/",sampleID,".rds"))


function_image_fixer(Spatial_Data) -> Spatial_Data



Spatial_Data[["DWLS"]] <- dwls_data


DefaultAssay(Spatial_Data) <- "DWLS"

cell_types_all=rownames(Spatial_Data)

wp=seurat_plotting()

ggsave(paste0("results/",sampleID,"/deconvolution/dwls/",sampleID,"-",scrnaID,"-dwls.pdf"),wp,height=18,width=8)