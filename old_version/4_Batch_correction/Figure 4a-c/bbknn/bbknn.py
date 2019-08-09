## generate bbknn coordinates ##

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



samples=['HCC1395','HCC1395BL','Mix']
datasets=['10x_LLU','10x_NCI','C1_FDA','C1_LLU','WaferGen']
datasets2=['10x_Mix10_LLU','10x_Mix5_NCI','10x_Mix5_F_NCI','10x_Mix5_F2_NCI']
hvgs=[100,500,1000,2000,4000]
num_pcs=20


for sample in samples:
	input=[]
	holder=[]
	datasets_chosen=datasets
	if sample=='Mix':
	 datasets_chosen=datasets2
	for dataset in datasets_chosen:
	 input.append(‘../manuscript_code/4_Batch_correction/Figure_4a-c/scanorama/data/'+sample+'/'+dataset+'.txt.gz')
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
	for hvg in hvgs:
	 filter_results = sc.pp.filter_genes_dispersion(adata.X, n_top_genes=hvg)
	 adata_new = adata[:, filter_results.gene_subset]
	 sc.pp.log1p(adata_new)
	 sc.pp.scale(adata_new, max_value=10)
	 sc.tl.pca(adata_new)
	 adata_new.obsm['X_pca'] *= -1
	 adata_bbknn = bbknn.bbknn(adata_new, neighbors_within_batch=4, n_pcs=num_pcs, trim=50, save_knn=True, copy=True)
	 sc.tl.umap(adata_bbknn)
	 umap_coord=pd.DataFrame(adata_bbknn.obsm['X_umap'])
	 sample_id=pd.DataFrame(adata_bbknn.obs['sample'])
	 umap_coord.index=sample_id.index
	 umap_coord=pd.concat([umap_coord,sample_id],axis=1)
	 output=‘../manuscript_code/4_Batch_correction/Figure_4a-c/bbknn/'+sample+'/HVG_'+str(hvg)
	 if not os.path.exists(output):
	  os.makedirs(output)
	 umap_coord.to_csv(output+'/umap_coord_bbknn.csv')




