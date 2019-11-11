import numpy as np
from scanorama import *
from sklearn.preprocessing import normalize, LabelEncoder
import sys

from process import load_names, merge_datasets, save_datasets


NAMESPACE = 'all'

data_names = [
    'data/All/10x_LLU_A',
    'data/All/10x_NCI_A',
    'data/All/10x_NCI_M_A',
    'data/All/10x_LLU_B',
    'data/All/10x_NCI_B',
    'data/All/10x_NCI_M_B',
    'data/All/10x_Mix10_LLU',
    'data/All/10x_Mix5_NCI',
    'data/All/10x_Mix5_F_NCI',
    'data/All/10x_Mix5_NCI_M',
    'data/All/10x_Mix5_F_NCI_M',
    'data/All/10x_Mix5_F2_NCI_M',
    'data/All/C1_FDA_HT_A',
    'data/All/C1_LLU_A',
    'data/All/ICELL8_SE_A',
    'data/All/ICELL8_PE_A',
    'data/All/C1_FDA_HT_B',
    'data/All/C1_LLU_B',
    'data/All/ICELL8_SE_B',
    'data/All/ICELL8_PE_B',
]

if __name__ == '__main__':
    datasets, genes_list, n_cells = load_names(data_names)

    datasets_dimred, datasets, genes = correct(
        datasets, genes_list, ds_names=data_names, hvg=2000, 
        return_dimred=True
    )

    save_datasets(datasets, genes, data_names)
