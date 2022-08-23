

rule seuratdecon:
    input:
        scrna="scrna/{datafile}.rds",
        spatial="analyses/raw/{sample}.rds",
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
        workflow/scripts/sp-tangram-pdf.R --rds {input.rds} --csv {input.csv} --output {output} --sampleid {wildcards.sample}
        """   



rule tangram_gene:
    input:
        "scrna/{datafile}.h5ad",
        "data/{sample}/outs/filtered_feature_bc_matrix.h5"
    output:
        "results/{sample}/deconvolution/tangramgene/{sample}-{datafile}-tangramgene.pdf"
    threads: 5
    resources:
        mem_mb=get_mem_mb,
        gpu=1
    shell:
        """
        workflow/scripts/spatial-tangram-gene.py data/{wildcards.sample}/outs {input[0]} {output}
        """   