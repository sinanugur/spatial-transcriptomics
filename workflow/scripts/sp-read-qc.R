#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("--min.cells"),
    type = "integer", default = 3,
    help = "Min cells [default= %default]", metavar = "integer"
  ),
  optparse::make_option(c("--min.features"),
    type = "integer", default = 200,
    help = "Min features [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--data.dir"),
    type = "character", default = NULL,
    help = "Data directory", metavar = "character"
  ),
  optparse::make_option(c("--sampleid"),
    type = "character", default = NULL,
    help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--percent.mt"),
    type = "double", default = 10,
    help = "Mitochondria filtering percentage [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--mad.nCount"),
    type = "double", default = 3,
    help = "Filtering percentage [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--mad.nFeature"),
    type = "double", default = 3,
    help = "Filtering percentage [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--output"),
    type = "character", default = "output.rds",
    help = "Output RDS file name", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$data.dir) || is.null(opt$sampleid)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (data.dir and sampleid)", call. = FALSE)
}

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)


# nFeature_Spatial is the number of genes detected in each cell. nCount_Spatial is the total number of molecules detected within a cell.


scrna <- Load10X_Spatial(
  data.dir = opt$data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1"
)


# scrna <- CreateSeuratObject(counts = scrna.data, project = opt$sampleid, min.cells = opt$min.cells, min.features = opt$min.features)


scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
scrna[["percent.rp"]] <- PercentageFeatureSet(scrna, pattern = "^RP[SL]")

scrna <- subset(scrna, subset = nFeature_Spatial > 0 & nCount_Spatial > 0) # remove empty cells

output.dir <- paste0("results/", opt$sampleid, "/technicals/")
dir.create(output.dir, recursive = T)

VlnPlot(scrna, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.rp"), ncol = 4)


ggsave(paste0(output.dir, "before-qc-trimming-violinplot.pdf"), width = 10, height = 4)




# minCov=opt$minCov #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good.
# if(min(scrna$nCount_Spatial)>=minCov){
#   countLOW=min(scrna$nCount_Spatial)
# }else{
#   countLOW=quantile(scrna$nCount_Spatial, prob=0.01)
# }
# countHIGH=quantile(scrna$nCount_Spatial, prob=0.99)
# featureLOW=quantile(scrna$nFeature_Spatial, prob=0.01)


lower_bound_nCount_Spatial <- median(scrna$nCount_Spatial) - opt$mad.nCount * mad(scrna$nCount_Spatial, constant = 1)
upper_bound_nCount_Spatial <- median(scrna$nCount_Spatial) + opt$mad.nCount * mad(scrna$nCount_Spatial, constant = 1)


lower_bound_nFeature_Spatial <- median(scrna$nFeature_Spatial) - opt$mad.nFeature * mad(scrna$nFeature_Spatial, constant = 1)
upper_bound_nFeature_Spatial <- median(scrna$nFeature_Spatial) + opt$mad.nFeature * mad(scrna$nFeature_Spatial, constant = 1)



## subset
# scrna <- subset(scrna, subset = nFeature_Spatial > featureLOW & nCount_Spatial > countLOW  & nCount_Spatial < countHIGH & opt$percent.mt < opt$percent.mt)

scrna <- subset(scrna, subset = nFeature_Spatial > lower_bound_nFeature_Spatial & nFeature_Spatial < upper_bound_nFeature_Spatial & nCount_Spatial > lower_bound_nCount_Spatial & nCount_Spatial < upper_bound_nCount_Spatial & percent.mt < opt$percent.mt)

VlnPlot(scrna, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.rp"), ncol = 4)
ggsave(paste0(output.dir, "after-qc-trimming-violinplot.pdf"), width = 10, height = 4)


# output.dir=paste0("analyses/raw/")
# dir.create(output.dir,recursive = T)

scrna <- SCTransform(scrna, assay = "Spatial", verbose = FALSE)
DefaultAssay(scrna) <- "Spatial"

scrna@meta.data$orig.ident <- opt$sampleid

saveRDS(scrna, file = opt$output)