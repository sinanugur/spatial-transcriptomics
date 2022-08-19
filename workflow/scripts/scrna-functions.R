#!/usr/bin/env Rscript



function_pca_dimensions=function(Spatial_Data){
  
  pct <- Stdev(object = Spatial_Data, reduction = "pca") / sum(Stdev(object = Spatial_Data, reduction = "pca")) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
dimensionReduction <- min(co1, co2) 
  
}

function_gofunc=function(df,algorithm="weight01",statistics="ks",mapping="org.Hs.eg.db",ID="symbol",ontology = "BP") {
  

geneList <- df$p_val
names(geneList) <- df$gene
# Create topGOData object
GOdata <- new("topGOdata",
              ontology = ontology,
              allGenes = geneList,
              geneSelectionFun = function(x)x,
              annot = annFUN.org, mapping = mapping,ID=ID,nodeSize = 10)

resultsKS <- runTest(GOdata,algorithm = algorithm,statistic = statistics)


tab <- GenTable(GOdata, raw.p.value = resultsKS, topNodes = length(resultsKS@score), numChar = 120)

return(tab)
}



function_image_fixer=function(Spatial_Data,sampleID) {

  IMAGE=Read10X_Image(image.dir=paste0("data/",sampleID,"/outs/spatial"),image.name="tissue_fixed.png")


IMAGE@coordinates[Spatial_Data@images$slice1@coordinates %>% rownames(),] -> IMAGE@coordinates


Spatial_Data@images$"image" <- IMAGE

Spatial_Data@images$"image"@assay <- "Spatial"
Spatial_Data@images$"image"@assay <- "Spatial"

Spatial_Data@images$"image"@key <- paste0("image","_")

return(Spatial_Data)

}

seurat_plotting=function() {

wp=Seurat::SpatialFeaturePlot(
  object = Spatial_Data,
  features = cell_types_all,alpha = c(0.7, 1),pt.size.factor = 1.5,ncol=2,images=paste0("image")) & 
  scale_fill_viridis() &   
  theme(legend.title = element_text(size=4.5),legend.key.size = unit(0.5,"cm"),
  legend.text = element_text(size=3),legend.margin=margin(t = 0,b = 0.1, unit='cm'),plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

  return(wp)

}