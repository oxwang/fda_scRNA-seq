import numpy as np
from scanorama import *
from sklearn.preprocessing import normalize, LabelEncoder
import sys

from process import load_names, merge_datasets, save_datasets


NAMESPACE = 'tian'

data_names = [
    'data/Tian/sc_10x_3cl',
    'data/Tian/sc_10x_5cl',
    'data/Tian/sc_CEL-seq2',
    'data/Tian/sc_Drop-seq'
]

if __name__ == '__main__':
    datasets, genes_list, n_cells = load_names(data_names)

    datasets_dimred, datasets, genes = correct(
        datasets, genes_list, ds_names=data_names, hvg=2000, 
        return_dimred=True
    )

    save_datasets(datasets, genes, data_names)
