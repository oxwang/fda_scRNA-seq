python_path=/home/genomics/bioinformatics/app/miniconda/miniconda3/bin/python
$python_path

import pandas as pd
import numpy as np
import scanpy.api as sc
import anndata
import os
import subprocess
import bbknn
from scipy import sparse


## We used scenario 2 as an example ##

samples=['scenario_2']
datasets=['10x_LLU','10x_NCI','C1_FDA_HT','C1_LLU','iCell8_SE']
hvg=2000

for sample in samples:
	input=[]
	holder=[]
	datasets_chosen=datasets
	for dataset in datasets_chosen:
	 input.append('/.../4_Batch_correction/Data/'+sample+'/'+dataset+'.txt.gz')
	 counts=pd.read_table(input[-1],index_col=0)
	 holder.append(anndata.AnnData(X=counts.values).T)
	 holder[-1].X=sparse.csr_matrix(holder[-1].X)
	 holder[-1].var_names=counts.index
	 holder[-1].obs_names=counts.columns
	 holder[-1].obs['sample']=dataset
	 sc.pp.filter_cells(holder[-1], min_genes=200)
	 sc.pp.filter_genes(holder[-1], min_cells=3)
	adata = holder[0].concatenate(holder[1:], join='outer')
	adata.X = adata.X.tocsr()
	adata.raw = sc.pp.log1p(adata, copy=True)
	sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e5)
	filter_results = sc.pp.filter_genes_dispersion(adata.X, n_top_genes=hvg)
	adata_new = adata[:, filter_results.gene_subset]
	sc.pp.log1p(adata_new)
	sc.pp.scale(adata_new, max_value=10)
	sc.tl.pca(adata_new)
	adata_new.obsm['X_pca'] *= -1
	adata_bbknn = bbknn.bbknn(adata_new, neighbors_within_batch=4, n_pcs=30, trim=50, copy=True)
	sc.tl.umap(adata_bbknn)
	umap_coord=pd.DataFrame(adata_bbknn.obsm['X_umap'])
	sample_id=pd.DataFrame(adata_bbknn.obs['sample'])
	umap_coord.index=sample_id.index
	umap_coord=pd.concat([umap_coord,sample_id],axis=1)
	output='/.../4_Batch_correction/Results/'+sample+'/HVG_'+str(hvg)+'/bbknn'
	if not os.path.exists(output):
	 os.makedirs(output)
	umap_coord.to_csv(output+'/umap_coord_bbknn.csv')

