
rule seurat_integrate:
    input:
        expand("analyses/raw/{sample}.rds",sample=files)
    output:
        rds="analyses/integration/seurat/" + integration_id + ".resolution-{res}.seurat.rds",
        #umap="results/integration/seurat/" + integration_id + ".umap_before_integration.pdf",
        umap="results/integration/seurat/" + integration_id + ".resolution-{res}.umap_after_integration.pdf"
    shell:
        """
        workflow/scripts/sp-seurat-integration.R --rds "{input}" --sampleid {integration_id} --umap.plot {output.umap} --output {output.rds} --resolution {wildcards.res}
        """

rule clustermarkers_seurat_integration:
    input:
        rds="analyses/integration/seurat/" + integration_id + ".resolution-{res}.seurat.rds"
    output:
        positive="results/integration/seurat/" + integration_id + ".resolution-{res}.positive-markers-forAllClusters.xlsx",
        allmarkers="results/integration/seurat/" + integration_id + ".resolution-{res}.all-markers-forAllClusters.xlsx"
    shell:
        """
        workflow/scripts/sp-find-markers.R --rds {input.rds} --resolution {wildcards.res} --sampleid {integration_id} --logfc.threshold {logfc_threshold} --test.use {test_use} --output.xlsx.positive {output.positive} --output.xlsx.all {output.allmarkers}
        """
