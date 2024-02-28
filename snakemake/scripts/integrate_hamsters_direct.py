import os
from pathlib import Path

import scvi
import GPUtil
import scanpy as sc
import matplotlib.pyplot as pl
import pandas as pd
import numpy as np
from matplotlib.cm import get_cmap

### CHECK GPU AVAILABILITY ###
import torch
print('GPU available:', GPUtil.getAvailable())
print('CUDA available:', torch.cuda.is_available())
print('CUDA device:', torch.cuda.get_device_name(0))

### PARAMETERS ###
subset = 'tcells'
max_epochs = 100000
n_latent = 30
n_layers = 2
for subset in ['endothelial', 'neutrophils', 'macrophages', 'tcells']:
    print(f'Running for {subset}')
    pr_timecourse = f'/g/huber/users/peidli/data/PanCov19/PhoRob_timecourse_{subset}.h5'
    ma_timecourse = f'/g/huber/users/peidli/data/PanCov19/MesAur_timecourse_{subset}.h5'
    output = f'/g/huber/users/peidli/data/PanCov19/Cov19Hamster_timecourse_integrated_{subset}.h5'


    ### PREPARATION ###
    # Load data
    print('Loading data...')
    adata = sc.read_h5ad(pr_timecourse)
    adata.uns['name'] = 'Dwarfhamster'
    adata.obs['organism'] = 'Dwarfhamster'
    adata.obs['dosage'] = ['high dose' if 'hd' in x else 'low dose' if 'ld' in x else 'no dose' for x in adata.obs['orig.ident']]
    adata.obs['time'] = [x.split('_')[-1] for x in adata.obs.timepoint]
    adata.obs['hamster'] = [f'{x}_pr' for x in adata.obs.hamster]

    bdata = sc.read_h5ad(ma_timecourse)
    bdata.uns['name'] = 'Goldhamster'
    bdata.obs['organism'] = 'Goldhamster'
    bdata.obs['dosage'] = ['no dose' if x=='d0' else 'high dose' for x in bdata.obs['timepoint']]  # high dose
    bdata.obs['time'] = [x.upper() for x in bdata.obs.timepoint]
    bdata.obs['timepoint'] = [f'hd_{x}' if x!='D0' else x for x in bdata.obs.time]
    bdata.obs['hamster'] = [f'{x}_ma' for x in bdata.obs.hamster]

    # concat data and define features
    print('Concating data...')
    adatas = [adata, bdata]
    superdata = sc.concat({adata.uns['name']: adata for adata in adatas}, label='batch', index_unique='-', join='outer')  # outer join to keep all genes
    sc.pp.highly_variable_genes(superdata, flavor='seurat_v3', subset=False, n_top_genes=2000, batch_key='batch', layer='counts')
    subdata = superdata[:, superdata.var.highly_variable].copy()

    ### INTEGRATION ###
    # Prepare data for scVI, define covariates to remove
    print(f'Running scVI for {subset} ...')
    scvi.model.SCVI.setup_anndata(subdata, layer='counts', batch_key="batch",
                            categorical_covariate_keys=None,
                            continuous_covariate_keys=None)

    # Setup and train model (if needed)
    model = scvi.model.SCVI(subdata, n_layers=n_layers, n_latent=n_latent)
    model.train(max_epochs=max_epochs,
                train_size=0.9, 
                validation_size=0.05,
                early_stopping=True)
    print('Finished training')

    # save resulting model
    model.save(f'../models/Gold_Dwarf_timecourse_scvi_snake_{subset}', overwrite=True)

    # plot loss
    pl.plot(model.history['train_loss_epoch'])
    pl.xlabel('epochs')
    pl.ylabel('training loss')
    pl.grid(True)
    pl.savefig(f'../models/Gold_Dwarf_timecourse_scvi_snake_{subset}/loss_plot.pdf', bbox_inches='tight')
    pl.close()

    # build scVI embedding
    latent = model.get_latent_representation(subdata)
    superdata.obsm["X_scVI"] = latent
    sc.pp.neighbors(superdata, use_rep="X_scVI", key_added='neighbors_scVI')
    sc.tl.leiden(superdata, neighbors_key='neighbors_scVI', key_added='leiden_scVI')
    sc.tl.umap(superdata, min_dist=0.3, neighbors_key='neighbors_scVI')
    superdata.obsm['X_umap_scVI'] = superdata.obsm['X_umap']

    ### OTHER EMBEDDINGS ###

    # harmony
    print('Running harmony...')
    sc.pp.pca(subdata)  # always run pca before you do harmony
    sc.external.pp.harmony_integrate(subdata, key='batch')
    superdata.obsm['X_pca_harmony'] = subdata.obsm['X_pca_harmony']
    sc.pp.neighbors(superdata, use_rep="X_pca_harmony", key_added='neighbors_harmony')
    sc.tl.leiden(superdata, neighbors_key='neighbors_harmony', key_added='leiden_harmony')
    sc.tl.umap(superdata, min_dist=0.3, neighbors_key='neighbors_harmony')
    superdata.obsm['X_umap_harmony'] = superdata.obsm['X_umap']

    # without integration
    print('Raw embedding...')
    sc.pp.pca(superdata, use_highly_variable=True)
    sc.pp.neighbors(superdata)
    sc.tl.umap(superdata)
    superdata.obsm['X_umap_no_integration'] = superdata.obsm['X_umap']

    # # set colors
    # print('Setting colors...')
    # superdata.uns['dosage_colors'] = ['tab:red', 'tab:blue', 'tab:grey']  # high, low, no dose
    # superdata.uns['timepoint_colors'] = ['tab:grey', *get_cmap('Reds')(np.linspace(0.2,1,4)), *get_cmap('Blues')(np.linspace(0.2,1,4))]  # D0, hd D2/3/5/E14, ld D2/D3
    # superdata.uns['time_colors'] = get_cmap('viridis')(np.linspace(0,1,5)) # D0, D2, D3, D5, e14
    # superdata.uns['batch_colors'] = ['tab:green', 'tab:orange']
    # superdata.uns['organism_colors'] = ['tab:green', 'tab:orange']
    # hamster_colors = [None]*(6)
    # hamster_colors[::2] = get_cmap('Oranges')(np.linspace(0.2,0.8,3))
    # hamster_colors[1::2] = get_cmap('Greens')(np.linspace(0.2,0.8,3))
    # superdata.uns['hamster_colors'] = hamster_colors

    # save integrated superdata
    print('Writing superdata object...')
    superdata.write(output)
