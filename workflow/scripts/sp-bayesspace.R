#!/usr/bin/env Rscript

option_list = list(
  optparse::make_option(c("--pca.dimension"), type="integer", default=15, 
              help="PCA dimensions [default= %default]", metavar="integer"),
    optparse::make_option(c("--data.dir"), type="character", default=NULL, 
              help="Visium data directory", metavar="character"),
    optparse::make_option(c("--output.sce"), type="character", default="output.sce.rds", 
              help="Output RDS sce file name", metavar="character"),
    optparse::make_option(c("--output.enhanced"), type="character", default="output.enhanced.rds", 
              help="Output RDS enhanced file name", metavar="character"),
    optparse::make_option(c("--qplot"), type="character", default="qplot.pdf", 
              help="Qplot filename", metavar="character"),
    optparse::make_option(c("--cluster.plot"), type="character", default="cluster.plot.pdf", 
              help="Cluster plot filename", metavar="character"),
    optparse::make_option(c("--n.cluster"), type="character", default=NULL, 
              help="Number of clusters, if not given, autoselect", metavar="integer"),
    optparse::make_option(c("--sampleid"), type="character", default=NULL, 
              help="Sample ID", metavar="character")
    


)
 
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$data.dir)){
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


sce <- readVisium(opt$data.dir)


set.seed(102)
sce <- spatialPreprocess(sce, platform="Visium", 
                              n.PCs=25, n.HVGs=2000, log.normalize=TRUE)

sce <- qTune(sce, qs=seq(2, 15), platform="Visium", d=opt$pca.dimension,gamma=3)
qPlot(sce)
ggsave(filename=opt$qplot)

n=opt$n.cluster
if(is.null(n)) {

test=function(x,y) { angles=atan2(y,x)*180/pi
  
  return(angles[1] - angles[2])
  }

n=attr(sce,"q.logliks") %>% as.data.frame() %>% select(q,loglik) %>% mutate(lo=(-1*loglik)/1000) %>% mutate(x1=q,y1=lo,x2=lag(q,order_by = q),y2=lag(lo,order_by = q),x3=lead(q,order_by = q),y3=lead(lo,order_by = q)) %>% 
rowwise() %>% mutate(ang=test(c(x2-x1,x3-x1),c(y2-y1,y3-y1))) %>% ungroup() %>% dplyr::filter(!is.na(ang)) %>%
{if(any(.$ang < 0)) dplyr::filter(.,row_number() < dplyr::first(which(ang < 0))) else .} %>% dplyr::filter(ang >0) %>% dplyr::slice_min(order_by = ang,n=1) %>% pull(q)


if(length(n) == 0) { n=7}

}


print("cluster:")
print(n)
set.seed(149)
sce <- spatialCluster(sce, q=n, platform="Visium", d=opt$pca.dimension,
                           init.method="mclust", model="t", gamma=3,
                           nrep=10000, burn.in=100,
                           save.chain=TRUE)


sce$sample_name <- opt$sampleid


set.seed(149)
sce_enhanced <- spatialEnhance(sce, q=n, platform="Visium", d=opt$pca.dimension,
                                    model="t", gamma=3,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=10000, burn.in=100,
                                    save.chain=TRUE)


wrap_plots(clusterPlot(sce)  + scale_x_reverse() + scale_y_reverse(),
clusterPlot(sce_enhanced,is.enhanced = T) + scale_x_reverse() + scale_y_reverse()) -> wp

ggsave(filename=opt$cluster.plot,wp)


sce_enhanced <- enhanceFeatures(sce_enhanced, sce,
                                     feature_names=domanska_markers,model = "xgboost",
                                     nrounds=0)


sce_enhanced$sample_name <- opt$sampleid

saveRDS(sce,opt$output.sce)
saveRDS(sce_enhanced,opt$output.enhanced)

