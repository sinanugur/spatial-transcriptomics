rule bayesspace:
    input:
        data_directory + "/{sample}/outs/filtered_feature_bc_matrix.h5"
    output:
        sce="analyses/bayesspace/{sample}.sce.rds",
        #enhanced="analyses/bayesspace/{sample}.sce.enhanced.rds",
        clusterplot="results/{sample}/bayesspace/{sample}.cluster_plot.pdf",
        qplot="results/{sample}/bayesspace/{sample}.qplot.pdf"
    params:
        d=data_directory + "/{sample}/outs/"
    shell:
        "workflow/scripts/sp-bayesspace.R  --sampleid {wildcards.sample} --pca.dimension {bayesspace_pca_dimension} --n.cluster {bayesspace_n_clusters} --data.dir {params.d} --output.sce {output.sce}  --cluster.plot {output.clusterplot} --qplot {output.qplot}"


rule bayesspace_enhance:
    input:
        sce="analyses/bayesspace/{sample}.sce.rds",
    output:
        enhanced="analyses/bayesspace/{sample}.sce.enhanced.rds"
    params:
        d=data_directory + "/{sample}/outs/"
    shell:
        "workflow/scripts/sp-bayesspace.R  --sampleid {wildcards.sample} --pca.dimension {bayesspace_pca_dimension} --n.cluster {bayesspace_n_clusters} --data.dir {params.d} --output.sce {output.sce}  --cluster.plot {output.clusterplot} --qplot {output.qplot}"



rule bayesspaceplots:
    input:
        "analyses/bayesspace/{sample}.sce.enhanced.rds"
    output:
        directory("results/{sample}/bayesspace/plots/")
    shell:
        "workflow/scripts/sp-bayesspace-feature-plot.R --input {input} --sampleid {wildcards.sample}"


