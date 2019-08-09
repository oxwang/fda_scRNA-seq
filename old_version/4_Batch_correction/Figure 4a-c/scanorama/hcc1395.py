import numpy as np
from sklearn.preprocessing import normalize, LabelEncoder
import sys

from process import load_names, merge_datasets, save_datasets
from scanorama import correct, visualize, process_data
from scanorama import dimensionality_reduce

NAMESPACE = 'hcc1395'

data_names = [
    'data/HCC1395/10x_LLU',
    'data/HCC1395/10x_NCI',
    'data/HCC1395/C1_FDA',
    'data/HCC1395/C1_LLU',
    'data/HCC1395/WaferGen',
]

if __name__ == '__main__':
    datasets, genes_list, n_cells = load_names(data_names)
    datasets, genes = correct(datasets, genes_list)
    datasets = [ normalize(ds, axis=1) for ds in datasets ]
    datasets_dimred = dimensionality_reduce(datasets)

    save_datasets(datasets, genes, data_names)
    save_datasets(datasets_dimred, genes, data_names)
