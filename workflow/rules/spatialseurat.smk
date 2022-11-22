from collections import defaultdict
from yaml import load
import os




rule rds:
    input:
        #directory(data_directory + "/{sample}/outs/")
        data_directory + "/{sample}/outs/filtered_feature_bc_matrix.h5"
    output:
        "analyses/raw/{sample}.rds"
    params:
        d=data_directory + "/{sample}/outs/"
    shell:
        "workflow/scripts/sp-read-qc.R --data.dir {params.d} --output {output} --sampleid {wildcards.sample} --percent.mt {percent_mt} --mad.nFeature {mad_nFeature} --mad.nCount {mad_nCount}"



rule imagefix:
    input:
        data_directory + "/{sample}/outs/spatial/tissue_hires_image.png"
    output:
        data_directory + "/{sample}/outs/spatial/tissue_fixed.png"

    shell:
        """
        #convert {input} -colorspace HCL -channel R -evaluate set 67% +channel -colorspace sRGB {output}
        workflow/scripts/saturation 0.4 {input} {output}
        """

rule spatialfeatureplot:
    input:
        rds="analyses/raw/{sample}.rds"
    output:
        heatmap="results/{sample}/technicals/SpatialFeature_QC.pdf",
    shell:
        """
        workflow/scripts/sp-spatialfeatureplot.R --rds {input.rds} --output {output}

        """

rule imagetissue:
    input:
        data_directory + "/{sample}/outs/spatial/tissue_hires_image.png"
    output:
        "results/{sample}/TissueImage/{sample}.filtered.png",
        "results/{sample}/TissueImage/{sample}.original.png"

    shell:
        """
        #convert {input} -colorspace HCL -channel R -evaluate set 67% +channel -colorspace sRGB {output}
        workflow/scripts/saturation 0.4 {input} {output[0]}
        cp {input} {output[1]}
        """

rule clustree:
    input:
        "analyses/raw/{sample}.rds"
    output:
        clustree="results/{sample}/clusteringTree/clusteringTree-{sample}.pdf",
        heatmap="results/{sample}/technicals/DimHeatMap_plot.pdf",
        hvfplot="results/{sample}/technicals/highly-variable-features.pdf"
    shell:
        "workflow/scripts/sp-clusteringtree.R --rds {input} --output {output.clustree} --heatmap {output.heatmap} --hvfplot {output.hvfplot}"



rule normalization_pca_rds:
    input:
        "analyses/raw/{sample}.rds"
    output:
        "analyses/processed/{res}/{sample}.rds",
        "results/{sample}/resolution-{res}/{sample}.number-of-cells-per-cluster.xlsx"
    shell:
        "workflow/scripts/sp-normalization-pca.R --rds {input} --sampleid {wildcards.sample} --nfeature {highly_variable_features} --resolution {wildcards.res}"


rule umap_plot:
    input:
        rds="analyses/processed/{res}/{sample}.rds",
        imagefile=data_directory + "/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/resolution-{res}/{sample}.umap.pdf"
    shell:
        "workflow/scripts/sp-umap.R --rds {input.rds} --sampleid {wildcards.sample} --resolution {wildcards.res}"

rule clustermarkers:
    input:
        rds="analyses/processed/{res}/{sample}.rds"
    output:
        positive= "results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx",
        allmarkers= "results/{sample}/resolution-{res}/{sample}.all-markers-forAllClusters.xlsx"
    shell:
        """
        workflow/scripts/sp-find-markers.R --rds {input.rds} --resolution {wildcards.res} --logfc.threshold {logfc_threshold} --test.use {test_use} --output.xlsx.positive {output.positive} --output.xlsx.all {output.allmarkers}

        """


rule selected_markers_plots:
    input:
        rds="analyses/processed/{res}/{sample}.rds",
        imagefile=data_directory + "/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/resolution-{res}/selected-markers/selected-markers-dotplot.pdf",
        directory("results/{sample}/resolution-{res}/selected-markers/plots/")
    shell:
        "workflow/scripts/sp-selected-marker-plots.R --rds {input.rds} --resolution {wildcards.res} --sampleid {wildcards.sample}"


rule positive_markers_plots:
    input:
        rds="analyses/processed/{res}/{sample}.rds",
        excel="results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx",
        imagefile=data_directory + "/{sample}/outs/spatial/tissue_fixed.png"
    output:
        directory("results/{sample}/resolution-{res}/markers/")
    shell:
        "workflow/scripts/sp-marker-plots.R --rds {input.rds} --resolution {wildcards.res} --sampleid {wildcards.sample} --xlsx {input.excel}"


rule spatialfeatures:
    input:
        rds="analyses/raw/{sample}.rds",
        imagefile= data_directory + "/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/spatial-markers/{sample}.spatial_markers.xlsx",
        directory("results/{sample}/spatial-markers/plots")
    shell:
        """
        mkdir -p {output[1]}
        workflow/scripts/sp-spatial-markers.R --rds {input.rds} --sampleid {wildcards.sample}  --output {output[0]} --selection.method {selection_method}

        """

rule h5ad:
    input:
        rds="analyses/processed/{res}/{sample}.rds"
    output:
        "analyses/h5ad/{res}/{sample}.h5ad"
    shell:
        "workflow/scripts/sp-convert-to-h5ad.R --rds {input.rds} --output {output} --type Spatial"


rule celltype:
    input:
        "analyses/h5ad/{res}/{sample}.h5ad"
    
    output:
        outputdir=directory("analyses/celltypist/{res}/{sample}/"+ f"{celltypist_model}"),
        predicted="analyses/celltypist/{res}/{sample}/" + f"{celltypist_model}" + "/predicted_labels.csv",
        dotplot="results" + "/{sample}/resolution-{res}" +  "/celltype_annotation/annotation." + f"{celltypist_model}" +  ".dotplot.pdf",
        xlsx="results" + "/{sample}/resolution-{res}" + "/celltype_annotation/cluster_annotation_table." + f"{celltypist_model}" + ".xlsx"
        
    shell:
        """
        workflow/scripts/sp-celltypist.py {input} {output.dotplot} {output.outputdir} {output.xlsx} {celltypist_model}
        """