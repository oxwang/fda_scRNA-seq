import numpy as np
from scanorama import *
from sklearn.preprocessing import normalize, LabelEncoder
import sys

from process import load_names, merge_datasets, save_datasets


NAMESPACE = 'hcc1395'

data_names = [
    'data/HCC1395/10x_LLU',
    'data/HCC1395/10x_NCI',
    'data/HCC1395/C1_FDA_HT',
    'data/HCC1395/C1_LLU',
    'data/HCC1395/iCell8_SE',
]

if __name__ == '__main__':
    datasets, genes_list, n_cells = load_names(data_names)

    datasets_dimred, datasets, genes = correct(
        datasets, genes_list, ds_names=data_names, hvg=2000, 
        return_dimred=True
    )

    save_datasets(datasets, genes, data_names)
