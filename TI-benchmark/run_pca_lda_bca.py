#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run a PCA, LDA, and BCA on input data. Assumes input data is cells by genes.

Usage:
  run_pca_lda_bca.py <input_file>

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from docopt import docopt

from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA

f = "data/dynverse_data_v2/real/gold/hematopoiesis-gates_olsson.expression"

def run_PCA(df):
    num_dim = 20
    pcs = PCA(num_dim,iterated_power=10).fit_transform(df.values)
    return pd.DataFrame(pcs, index=df.index.values, columns=['PC'+str(i+1) for i in range(pcs.shape[1])])

def PCA_processed_LDA(df, label_df):
    num_dim = 20 # same as in structdr paper
    pcs = PCA(num_dim,iterated_power=10).fit_transform(df.values)
    Y_lda = LDA().fit_transform(pcs, dataanno['group_id'])
    return pd.DataFrame(Y_lda, index=df.index.values, columns=['LD'+str(i+1) for i in range(Y_lda.shape[1])])

def BCA(df, label_df):
    df = df - df.mean() # mean center genes of data
    X = df.to_numpy()
    df['cell_type'] = label_df['group_id'].to_list()
    M = df.groupby('cell_type').mean().to_numpy()
    D = np.diag(df['cell_type'].value_counts().sort_index().to_numpy())
    U, S, V_T = np.linalg.svd(D @ M, full_matrices=False)
    Y =  X @ V_T.T  # assumes X is cells by feats
    return pd.DataFrame(Y, index=df.index.values, columns=['BCC'+str(i+1) for i in range(Y.shape[1])])


# python ../run_structdr.py --input ./data/dynverse_data_v2/synthetic/splatter/tree_1.expression.dim20_k10_reg500_n0t12_pca.       graphdr --output ./data/dynverse_data_v2/synthetic/splatter/tree_1.expression.outdir/dim20_k10_reg500_n0t12_pca.graphdr  --ndim  6  --n_jobs 1 --anno_file ./data/dynverse_data_v2/synthetic/splatter/tree_1.anno --extract_backbone

if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.01')
    f = arguments['<input_file>']
    data = pd.read_csv(f,sep='\t',index_col=0)
    dataanno = pd.read_csv(f.replace('expression','anno'),sep='\t',index_col=0)

    Y_pca = run_PCA(data)
    Y_pca.to_csv(f+'.pca20',index_label=False,sep='\t')

    Y_lda = PCA_processed_LDA(data, dataanno)
    Y_lda.to_csv(f+f".lda{Y_lda.shape[1]}", index_label=False,sep='\t')

    Y_bca = BCA(data, dataanno).iloc[:, :-1] # for some reason SVD does not return num_classes - 1, but num_classes
    Y_bca.to_csv(f+f".bca{Y_bca.shape[1]}", index_label=False,sep='\t')
