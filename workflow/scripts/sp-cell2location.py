#!/usr/bin/env python

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi
import sys

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
from cell2location.models import RegressionModel
import torch



print(sys.argv[1])
adata_vis = sc.read(sys.argv[1])
adata_ref = sc.read(sys.argv[2])

adata_ref.X=adata_ref.raw.X.copy()



try:
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                            # 10X reaction / sample / batch
                            batch_key='orig.ident',
                            # cell type, covariate used for constructing signatures
                            labels_key='seurat_clusters'
                        )
except:
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                            # 10X reaction / sample / batch
                            #batch_key='orig.ident',
                            # cell type, covariate used for constructing signatures
                            labels_key='seurat_clusters'
                        )



# create the regression model
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
#mod.view_anndata_setup()



mod.train(max_epochs=250, use_gpu=torch.cuda.is_available())


adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': torch.cuda.is_available()}
)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="orig.ident")

mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
#mod.view_anndata_setup()

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=torch.cuda.is_available())

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(10)
#plt.legend(labels=['full data training']);

adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': torch.cuda.is_available()}
)

df=adata_vis.obsm['q05_cell_abundance_w_sf']

df.columns = df.columns.str.replace("q05cell_abundance_w_sf_","")

df.to_csv(sys.argv[3])