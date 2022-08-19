

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

