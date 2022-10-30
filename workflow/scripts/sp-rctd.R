#!/usr/bin/env Rscript
option_list = list(


    optparse::make_option(c("--sprds"), type="character", default=NULL, 
              help="Processed rds file of a spatial object", metavar="character"),
    optparse::make_option(c("--scrds"), type="character", default=NULL, 
              help="Processed rds file of a single cell object", metavar="character"),
        optparse::make_option(c("--output"), type="character", default="output.csv", 
              help="Output CSV file name", metavar="character")



)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$sprds)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (sprds)", call.=FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(spacexr)


source("workflow/scripts/scrna-functions.R")

Spatial_Data=readRDS(opt$sprds)
scrna_data=readRDS(opt$scrds)


cell_types_all=Idents(scrna_data) %>% unique() %>% as.character() %>% make.names()



counts <- scrna_data@assays$RNA@counts
colnames(counts) <- colnames(scrna_data)
meta_data <- data.frame(scrna_data@meta.data)
cell_types <- meta_data$seurat_clusters %>% make.names()
names(cell_types) <- rownames(scrna_data@meta.data)
cell_types <- as.factor(cell_types)
nUMI_df <- data.frame(colSums(scrna_data@assays$RNA@counts))
nUMI <- nUMI_df$colSums.scrna_data.assays.RNA
names(nUMI) <- rownames(nUMI_df)



reference <- Reference(counts, cell_types, nUMI)




coords <- data.frame(colnames(Spatial_Data))
colnames(coords) <- 'barcodes'
coords$xcoord <- seq_along(colnames(Spatial_Data))
coords$ycoord <- seq_along(colnames(Spatial_Data))
counts <- data.frame(Spatial_Data@assays$Spatial@counts) # load in counts matrix
colnames(counts) <- colnames(Spatial_Data)
coords <- data.frame(colnames(Spatial_Data))
colnames(coords) <- 'barcodes'
coords$xcoord <- seq_along(colnames(Spatial_Data))
coords$ycoord <- seq_along(colnames(Spatial_Data))
rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

puck<- SpatialRNA(coords = coords,counts = counts,nUMI = nUMI)

myRCTD <- create.RCTD(puck,reference,max_cores = 5,UMI_min = 20,CELL_MIN_INSTANCE=5)
myRCTD <- run.RCTD(myRCTD,doublet_mode = "full")


Spatial_Data@meta.data <- Spatial_Data@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(myRCTD@results$weights %>% as.data.frame() %>% rownames_to_column("barcodes"), by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")


#saveRDS(Spatial_Data,file = paste0("rds_rctd/",scrnaID,"/",sampleID,".rds"))


write.csv(myRCTD@results$weights %>% as.data.frame(),file=opt$output)

