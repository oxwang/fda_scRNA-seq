import numpy as np
from sklearn.preprocessing import normalize, LabelEncoder
import sys

from process import load_names, merge_datasets, save_datasets
from scanorama import correct, visualize, process_data
from scanorama import dimensionality_reduce

data_names = [
    'data/Mix/10x_Mix10_LLU',
    'data/Mix/10x_Mix5_NCI',
    'data/Mix/10x_Mix5_F_NCI',
    'data/Mix/10x_Mix5_F2_NCI',
]

if __name__ == '__main__':
    datasets, genes_list, n_cells = load_names(data_names)
    datasets, genes = correct(datasets, genes_list)
    datasets = [ normalize(ds, axis=1) for ds in datasets ]
    datasets_dimred = dimensionality_reduce(datasets)

    save_datasets(datasets, genes, data_names)
    save_datasets(datasets_dimred, genes, data_names)
