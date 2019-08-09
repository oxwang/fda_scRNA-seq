## generate umap coordinates ##

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



method=['unc','cca','mnn','limma','combat']
hvg=['HVG_100','HVG_500','HVG_1000','HVG_2000','HVG_4000']
cell=['HCC1395',"HCC1395BL","Mix"]

for i in range(3):
 for j in range(5):
  for k in range(5):
   input='/genomics/1_Projects/FDA_QC/Results/batch/Eva_1/Multi_data/'+cell[i]+'/'+'/'+hvg[j]+'/'+method[k]+'/'+method[k]+'_correct.csv'
   output='/genomics/1_Projects/FDA_QC/Results/batch/Eva_1/Multi_data/'+cell[i]+'/'+'/'+hvg[j]+'/'+method[k]+'/'+method[k]+'_umap.csv'
   counts=pd.read_table(input,index_col=0,sep=',')
   adata=anndata.AnnData(X=counts.values)
   if method[k]=='cca':
    adata.obsm['X_pca']=counts.values
   else:
    sc.pp.scale(adata,max_value=10)
    sc.tl.pca(adata)
    adata.obsm['X_pca']*=-1
   sc.pp.neighbors(adata,n_pcs=20,n_neighbors=20)
   sc.tl.umap(adata)
   umap_coord=pd.DataFrame(adata.obsm['X_umap'])
   umap_coord.index=counts.index
   umap_coord.to_csv(output)













