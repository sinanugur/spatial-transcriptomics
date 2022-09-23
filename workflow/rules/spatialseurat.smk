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



rule bayesspace:
    input:
        data_directory + "/{sample}/outs/filtered_feature_bc_matrix.h5"
    output:
        sce="analyses/bayesspace/{sample}.sce.rds",
        enhanced="analyses/bayesspace/{sample}.sce.enhanced.rds",
        clusterplot="results/{sample}/bayesspace/{sample}.cluster_plot.pdf",
        qplot="results/{sample}/bayesspace/{sample}.qplot.pdf"
    params:
        d=data_directory + "/{sample}/outs/"
    shell:
        "workflow/scripts/sp-bayesspace.R --pca.dimension {bayesspace_pca_dimension} --data.dir {params.d} --output.sce {output.sce} --output.enhanced {output.enhanced} --cluster.plot {output.clusterplot} --qplot {output.qplot}"

rule bayesspaceplots:
    input:
        "analyses/bayesspace/{sample}.sce.enhanced.rds"
    output:
        directory("results/{sample}/bayesspace/plots/")
    shell:
        "workflow/scripts/sp-bayesspace-feature-plot.R --input {input} --sampleid {wildcards.sample}"

rule imagefix:
    input:
        data_directory + "/{sample}/outs/spatial/tissue_lowres_image.png"
    output:
        "{data_directory}{sample}/outs/spatial/tissue_fixed.png"

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
        rds="analyses/processed/{res}/{sample}.rds"
    output:
        "results/{sample}/resolution-{res}/{sample}.umap.pdf"
    shell:
        "workflow/scripts/sp-umap.R --rds {input.rds} --sampleid {wildcards.sample} --resolution {wildcards.res}"


    
rule clustermarkers:
    input:
        rds="analyses/processed/{res}/{sample}.rds",
        imagefile=data_directory + "/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx",
        "results/{sample}/resolution-{res}/{sample}.all-markers-forAllClusters.xlsx"
    shell:
        """
        workflow/scripts/sp-find-markers.R --rds {input.rds} --resolution {wildcards.res} --sampleid {wildcards.sample} --logfc.threshold {logfc_threshold} --test.use {test_use}
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
        excel="results/{sample}/resolution-{res}/{sample}.positive-markers-forAllClusters.xlsx"
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

