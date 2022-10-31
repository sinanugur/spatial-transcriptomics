


rule scttransform_spatial:
    input:
        spatial="analyses/raw/{sample}.rds",
    output:
        "analyses/sct/{sample}.rds"
    shell:
        "workflow/scripts/sp-seurat-sct.R --rds {input.spatial} --output {output} --assay Spatial"  


rule scttransform_scrna:
    input:
        scrna="scrna/{datafile}.rds",
    output:
        sct="analyses/sct_scrna/{datafile}.rds",
        xlsx="analyses/sct_scrna/{datafile}.cluster_markers.xlsx"
    shell:
        "workflow/scripts/sp-seurat-sct.R --rds {input.scrna} --output {output.sct} --assay RNA --output.xlsx {output.xlsx} "  

rule seuratdecon:
    input:
        scrna="analyses/sct_scrna/{datafile}.rds",
        spatial="analyses/sct/{sample}.rds",
        imagefile="data/{sample}/outs/spatial/tissue_fixed.png"
    output:
        "results/{sample}/deconvolution/seurat/{sample}-{datafile}-seurat.pdf"
    shell:
        """
        workflow/scripts/sp-seurat-decon.R --sprds {input.spatial} --scrds {input.scrna} --output {output}  --sampleid {wildcards.sample}
        """

rule spatial_to_h5ad:
    input:
        rds="analyses/raw/{sample}.rds"
    output:
        "analyses/h5ad/{sample}.h5ad"
    shell:
        "workflow/scripts/sp-convert-to-h5ad.R --rds {input.rds} --output analyses/h5ad/{wildcards.sample} --type Spatial"


rule scrna_to_h5ad:
    input:
        rds="scrna/{datafile}.rds"
    output:
        "scrna/{datafile}.h5ad"
    shell:
        "workflow/scripts/sp-convert-to-h5ad.R --rds {input.rds} --output scrna/{wildcards.datafile} --type RNA"



rule tangram_cluster:
    input:
        spatial="analyses/h5ad/{sample}.h5ad",
        scrna="scrna/{datafile}.h5ad"
    output:
        "analyses/tangram/{datafile}/{sample}.csv"
    threads: 3
    resources:
        mem_mb=get_mem_mb,
        gpu=1
    shell:
        """
        workflow/scripts/sp-tangram.py {input.spatial} {input.scrna} {output}
        """

rule tangram_pdf:
    input:
        rds="analyses/raw/{sample}.rds",
        csv="analyses/tangram/{datafile}/{sample}.csv"
    output:
        "results/{sample}/deconvolution/tangram/{sample}-{datafile}-tangram.pdf"
    shell:
        """
        workflow/scripts/sp-features-pdf.R --rds {input.rds} --csv {input.csv} --output {output} --sampleid {wildcards.sample}
        """   



rule tangram_gene:
    input:
        spatial="analyses/h5ad/{sample}.h5ad",
        scrna="scrna/{datafile}.h5ad"
    output:
        #"results/{sample}/deconvolution/tangramgene/{sample}-{datafile}-tangramgene.pdf"
        predicted="analyses/tangramgene/{datafile}/{sample}.predicted.csv",
        measured="analyses/tangramgene/{datafile}/{sample}.measured.csv",

    threads: 5
    resources:
        mem_mb=get_mem_mb,
        gpu=1
    shell:
        """
        workflow/scripts/sp-tangram-gene.py {input.spatial} {input.scrna} {output.predicted} {output.measured}
        """

rule tangram_gene_pdf:
    input:
        rds="analyses/raw/{sample}.rds",
        predicted="analyses/tangramgene/{datafile}/{sample}.predicted.csv",
        measured="analyses/tangramgene/{datafile}/{sample}.measured.csv"
    output:
        directory("results/{sample}/deconvolution/tangramgene/{datafile}/")
    shell:
        """
        workflow/scripts/sp-tangram-gene-pdf.R --rds {input.rds} --predicted {input.predicted} --measured {input.measured} --output.dir {output} --sampleid {wildcards.sample}
        """   

rule cell2location:
    input:
        spatial="analyses/h5ad/{sample}.h5ad",
        scrna="scrna/{datafile}.h5ad"
    output:
        "analyses/cell2location/{datafile}/{sample}.csv"
    threads: 5
    resources:
        mem_mb=get_mem_mb,
        gpu=1
    shell:
        """
        workflow/scripts/sp-cell2location.py {input.spatial} {input.scrna} {output}
        """   

rule cell2loc_pdf:
    input:
        rds="analyses/raw/{sample}.rds",
        csv="analyses/cell2location/{datafile}/{sample}.csv"
    output:
        "results/{sample}/deconvolution/cell2location/{sample}-{datafile}-cell2location.pdf"
    shell:
        """
        workflow/scripts/sp-features-pdf.R --rds {input.rds} --csv {input.csv} --output {output} --sampleid {wildcards.sample}
        """  

rule bayesspace_integrate:
    input:
        expand("analyses/bayesspace/{sample}.sce.rds",sample=files)
    output:
        rds="analyses/integration/bayesspace/" + integration_id + "_bayesspace.rds",
        umap="results/integration/bayesspace/" + integration_id + ".umap_before_integration.pdf",
        harmony="results/integration/bayesspace/" + integration_id + ".umap_after_integration.pdf",


    shell:
        """
        workflow/scripts/sp-bayesspace-integration.R --rds "{input}" --sampleid {integration_id} --output {output.rds} --umap.plot {output.umap} --harmony.plot {output.harmony}
        """


rule giotto_dwls:
    input:
        scrna="scrna/{datafile}.rds",
        spatial="analyses/raw/{sample}.rds"
    output:
        "analyses/dwls/{datafile}/{sample}.csv"
    shell:
        """
        workflow/scripts/sp-giotto-dwls.R --sprds {input.spatial} --scrds {input.scrna} --output {output}
        """


rule rctd_decon:
    input:
        scrna="scrna/{datafile}.rds",
        spatial="analyses/raw/{sample}.rds"
    output:
        "analyses/rctd/{datafile}/{sample}.csv"
    shell:
        """
        workflow/scripts/sp-rctd.R --sprds {input.spatial} --scrds {input.scrna} --output {output}
        """

rule rctd_pdf:
    input:
        rds="analyses/raw/{sample}.rds",
        csv="analyses/rctd/{datafile}/{sample}.csv"
    output:
        "results/{sample}/deconvolution/rctd/{sample}-{datafile}-rctd.pdf"
    shell:
        """
        workflow/scripts/sp-features-pdf.R --rds {input.rds} --csv {input.csv} --output {output} --sampleid {wildcards.sample}
        """  

rule dwls_pdf:
    input:
        rds="analyses/raw/{sample}.rds",
        csv="analyses/cell2location/{datafile}/{sample}.csv"
    output:
        "results/{sample}/deconvolution/cell2location/{sample}-{datafile}-cell2location.pdf"
    shell:
        """
        workflow/scripts/sp-features-pdf.R --rds {input.rds} --csv {input.csv} --output {output} --sampleid {wildcards.sample}
        """  