import numpy as np
from scanorama import *
from sklearn.preprocessing import normalize, LabelEncoder
import sys

from process import load_names, merge_datasets, save_datasets


NAMESPACE = 'mix'

data_names = [
    'data/Mix/10x_Mix10_LLU',
    'data/Mix/10x_Mix5_NCI',
    'data/Mix/10x_Mix5_F_NCI',
    'data/Mix/10x_Mix5_F2_NCI',
]

if __name__ == '__main__':
    datasets, genes_list, n_cells = load_names(data_names)

    datasets_dimred, datasets, genes = correct(
        datasets, genes_list, ds_names=data_names, hvg=2000, 
        return_dimred=True
    )

    save_datasets(datasets, genes, data_names)
