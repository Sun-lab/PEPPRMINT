
from __future__ import absolute_import, division, print_function, unicode_literals

import warnings

warnings.filterwarnings(action="ignore", category=DeprecationWarning)
warnings.filterwarnings(action="ignore", category=FutureWarning)
warnings.filterwarnings(action="ignore", category=PendingDeprecationWarning)

import os
import subprocess
import argparse
import itertools
import random
import math
import numpy as np
import pandas as pd
import sys

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Dense, Flatten, Reshape
from tensorflow.keras.layers import Dropout, concatenate

from sklearn.preprocessing import normalize
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score

from datetime import datetime
from scipy import stats
from collections import Counter

import matplotlib
matplotlib.use("Agg") # use a non-interactive backend
# matplotlib.use('macosx')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------
# utility function to encode peptides by one-hot coding
# -------------------------------------------------------------------------

AA_SYMOLS = ['A', 'R', 'N', 'D', 'C',
             'Q', 'E', 'G', 'H', 'I',
             'L', 'K', 'M', 'F', 'P',
             'S', 'T', 'W', 'Y', 'V', 'X']

def encode_peptide(peptide):
    m = list()
    for aa in peptide:
        channel = list()
        for symbol in AA_SYMOLS:
            if aa.upper() == symbol: channel.append(0.90)
            else: channel.append(0.05)
        m.append(channel)
    m = np.array(m)
    return m

# -------------------------------------------------------------------------
# main function 
# -------------------------------------------------------------------------

# input_test_pred = "riaz_peptide_mut.txt"
# input_alist     = "riaz_HLA_I_allele_list.txt"
# input_MHCseq    = "MHC_pseudo.txt"
# test_data_name  = "PEPPRMINT_Riaz_2017"
# model           = "PEPPRMINT_800_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split4_Aug5_iter10.h5"
# m_tag           = "800_split4"
# data_dir        = "../data/Riaz_2017/"
# info_dir        = "../data/NetMHCpan4_1_train/"
# fig_dir         = "../figures/PEPPRMINT"
# results_dir     = "../results/Riaz_2017"
# model_dir       = "../results/PEPPRMINT/"
# save_all_pred   = True
# input_pep_hla   = False
# encode9AA       = False
# neoantigen      = True

def main(input_test_pred, input_alist, 
    input_MHCseq, test_data_name,  model, m_tag, 
    data_dir, fig_dir, results_dir, info_dir, model_dir,
    save_all_pred = False, input_pep_hla = False, 
    encode9AA = False, neoantigen = False):

    saved_args = locals()
    print("saved_args is", saved_args)

    config = str(test_data_name)  + '_' +  str(m_tag) 
    model  = str(model) 

    # -----------------------------------------------------------------
    # read in data
    # -----------------------------------------------------------------

    print('read in data')
    file_test_pred = os.path.join(data_dir, input_test_pred)
    test_predX     = pd.read_csv(file_test_pred, sep='\t', header=0)

    print('data dimension:')
    print(test_predX.shape)

    print('test_predX[0:2,]:')
    print(test_predX.iloc[0:2,])

    # test_predX['binder'].value_counts()
    # test_predX['HLA'].value_counts()

    # -----------------------------------------------------------------
    # read in HLA allele information for each cell line
    # -----------------------------------------------------------------

    file_allele = os.path.join(input_alist)
    allele_list = pd.read_csv(file_allele, sep=' ', header=None)
    print(allele_list.shape)
    print(allele_list.iloc[0:3,])
    print(allele_list.iloc[1,1].split(","))

    allele_dic = {}

    for i in allele_list.index:
        allele_dic[allele_list.iloc[i,0]] = allele_list.iloc[i,1].split(",")

    first2pairs = {k: allele_dic[k] for k in list(allele_dic)[:2]}
    first2pairs

    # -----------------------------------------------------------------
    # read in HLA seqeunces
    # -----------------------------------------------------------------

    file_hla = os.path.join(info_dir, input_MHCseq)
    hla_seq  = pd.read_csv(file_hla, sep=' ', header=None)
    hla_seq.shape
    print(hla_seq.iloc[0:3,:])

    # -----------------------------------------------------------------
    # encode all the HLA seqeunces
    # -----------------------------------------------------------------

    hla_encoding = {}

    for i in hla_seq.index:
        hla_name = hla_seq.iloc[i,0]
        hla_aa   = hla_seq.iloc[i,1]
        hla_encoding[hla_name] = encode_peptide(hla_aa)

    print('finished hla_encode')

    # -----------------------------------------------------------------
    # For final evaluation of model: using input_test_pred
    # encode all possible HLA alleles for either binder or non-binder
    # -----------------------------------------------------------------

    test_idx     = list()
    test_hla     = list()
    test_encode  = list()
    test_encode_hla = list()

    if input_pep_hla:
        print('prepare pep and hla as separate input')

    for i in test_predX.index:
        pep_encode = encode_peptide(test_predX.iloc[i,0])
        cell_line  = test_predX.iloc[i,2]
        hla_allele = allele_dic[cell_line]
            
        for hla1 in hla_allele:
            hla_encode = hla_encoding[hla1]
            test_idx.append(i)
            test_hla.append(hla1)
                            
            if input_pep_hla:
                test_encode.append((pep_encode).flatten())
                test_encode_hla.append((hla_encode).flatten())
            else: 
                test_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())

    df_test = test_predX.loc[test_idx, :]
    col_nms = df_test.columns.values
    col_nms

    if neoantigen:
        col_nms[0:3] = ['peptide', 'key', 'sample']
    else:
        col_nms[0:3] = ['peptide', 'y_true', 'cell_line']
        
    df_test.columns = col_nms
    df_test

    test_encode      = np.asarray(test_encode)
    test_encode_hla  = np.asarray(test_encode_hla)

    # -----------------------------------------------------------------
    # load model for prediction
    # -----------------------------------------------------------------

    file_model = os.path.join(model_dir, model)
    model1 = tf.keras.models.load_model(file_model)

    # -----------------------------------------------------------------
    # Predict for weighted and unweighted models
    # -----------------------------------------------------------------
    if input_pep_hla: 
        test_pred = model1.predict([test_encode, test_encode_hla]) 
    else: 
        test_pred = model1.predict(test_encode) 

    df_test['y_pred'] = test_pred

    # -----------------------------------------------------------------
    # Take maximum for peptide/cell_line and calculate classification
    # -----------------------------------------------------------------

    df_test.iloc[0:10,]

    if encode9AA:
        var2grp = ['peptide_orig', 'y_true', 'cell_line']
    elif neoantigen: 
        var2grp = ['key', 'sample']
    else:
        var2grp = ['peptide', 'y_true', 'cell_line']

    df_max  = df_test.groupby(var2grp, as_index=False)['y_pred'].max()

    if not neoantigen:
        auc_roc = roc_auc_score(df_max["y_true"], df_max["y_pred"])
        pr_roc  = average_precision_score(df_max["y_true"], df_max["y_pred"])
        
        print('evaluation (of max if MA)')
        print('auc ROC, auc PR')
        print([round(auc_roc,3), round(pr_roc,3)])
        
        df_max.groupby('y_true').describe()
        
        print('')
        print('classificiation accuracy for cutoff 0.75:')
        print('')
        print(pd.crosstab(df_max['y_pred'] > 0.75, df_max['y_true']))
        
        print('')
        print('classificiation accuracy for cutoff 0.50:')
        print('')
        print(pd.crosstab(df_max['y_pred'] > 0.50, df_max['y_true']))
        print('')

    # -----------------------------------------------------------------
    # Save predictions
    # -----------------------------------------------------------------

    fnm_max = config +'_test_pred.txt'
    df_max.to_csv(os.path.join(results_dir, fnm_max), 
        sep='\t', index=False, float_format='%.3e')
    print('Done with predictions!')
        
    if save_all_pred:
        fnm  = config +'_test_pred_all.txt'
        df_test.to_csv(os.path.join(results_dir, fnm), 
            sep='\t', index=False, float_format='%.3e')
        subprocess.run(["gzip", os.path.join(results_dir, fnm)])

# -----------------------------------------------------------------
# parameters
# -----------------------------------------------------------------

parser = argparse.ArgumentParser(
    description='a pan-specific method for antigen presentation using deep learning.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "--input_test_pred", 
    type = str,
    dest = "input_test_pred",
    help = "input: file name of indpendent testing data. Delimited by a tab. "
)
parser.add_argument(
    "--input_alist", 
    type = str,
    dest = "input_alist",
    default = "../data/NetMHCpan4_1_train/allelelist.txt",
    help = "input: file name for list of HLA-I alleles for all the samples. Delimited by a single sapce."
)
parser.add_argument(
    "--input_MHCseq", 
    type = str,
    dest = "input_MHCseq",
    default = "MHC_pseudo.txt",
    help = "input: file with 34 pseudo-sequence of different HLA alleles. Delimited by a single space"
)
parser.add_argument(
    "--test_data_name", 
    type = str,
    dest = "test_data_name",
    help = "input: file name of testing data"
)
parser.add_argument(
    "--model", 
    type = str,
    dest = "model",
    help = "input: .h5 model which to make prediction with"
)
parser.add_argument(
    "--m_tag",
    type = str,
    dest = "m_tag",
    help = "short string identifying model which will be used to in file name of saved results"
)
parser.add_argument(
    "--data_dir", "-D",
    type = str,
    default = "../data/test_data",
    help = "path to directory where data are placed"
)
parser.add_argument(
    "--fig_dir", "-F",
    type = str,
    default = "../figures/test_data",
    help = "path to directory where figures are saved"
)
parser.add_argument(
    "--results_dir", "-R",
    type = str,
    default = "../results/test_data",
    help = "path to directory where results are saved"
)
parser.add_argument(
    "--info_dir", "-I",
    type = str,
    dest = "info_dir",
    default = "../data/NetMHCpan4_1_train",
    help = "directory where HLA allele informaton are placed"
)
parser.add_argument(
    "--model_dir",
    type = str,
    default = "../results/PEPPRMINT/",
    help = "path to directory where results are saved"
)
parser.add_argument(
    "--save_all_pred",
    dest = "save_all_pred",
    action = "store_true",
    default = False, 
    help = "option to save all predictions for all combinations of peptide and HLA in cell line"
)
parser.add_argument(
    "--input_pep_hla",
    dest = "input_pep_hla",
    action = "store_true",
    default = False, 
    help = "option to use v1 of model (i.e. peptide and hla are seperate inputs)"
)
parser.add_argument(
    "--encode9AA",
    dest = "encode9AA",
    action = "store_true",
    default = False, 
    help = "whether the testing data are encoded as 9AA sequence"
)
parser.add_argument(
    "--neoantigen",
    dest = "neoantigen",
    action = "store_true",
    default = False, 
    help = "whether the testing data for neoantigen prediction"
)


if __name__ == '__main__':
    arguments = parser.parse_args()
    main(**vars(arguments))

