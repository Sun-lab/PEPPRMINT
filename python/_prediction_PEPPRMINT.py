from __future__ import absolute_import, division, print_function, unicode_literals

import warnings

warnings.filterwarnings(action="ignore", category=DeprecationWarning)
warnings.filterwarnings(action="ignore", category=FutureWarning)
warnings.filterwarnings(action="ignore", category=PendingDeprecationWarning)

import os
import subprocess
import argparse
import numpy as np
import pandas as pd

import tensorflow as tf
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score

from datetime import datetime

import matplotlib

matplotlib.use("Agg")  # use a non-interactive backend
# matplotlib.use('macosx')

# -------------------------------------------------------------------------
# utility function to encode peptides by one-hot coding
# -------------------------------------------------------------------------

AA_SYMOLS = ['A', 'R', 'N', 'D', 'C',
             'Q', 'E', 'G', 'H', 'I',
             'L', 'K', 'M', 'F', 'P',
             'S', 'T', 'W', 'Y', 'V', 'X']


def encode_peptide(peptide):
    """ encode the peptide by one-hot encoding"""
    m = list()
    for aa in peptide:
        channel = list()
        for symbol in AA_SYMOLS:
            if aa.upper() == symbol:
                channel.append(0.90)
            else:
                channel.append(0.05)
        m.append(channel)
    m = np.array(m)
    return m


# -------------------------------------------------------------------------
# main function 
# -------------------------------------------------------------------------

# input_peptides  = "../data/Riaz_2017/riaz_peptide_mut.txt"
# input_alist     = "../data/Riaz_2017/riaz_HLA_I_allele_list.txt"
# input_mhc_seq   = "../data/NetMHCpan4_1_train/MHC_pseudo.dat"
# data_label      = "PEPPRMINT_Riaz_2017"
# model           = "MA_200_split0.h5"
# model_tag       = "200_split0"
# fig_dir         = "../figures/PEPPRMINT"
# results_dir     = "../results/Riaz_2017"
# save_all_pred   = True
# input_pep_hla   = False
# encode9AA       = False
# neoantigen      = True

def main(input_peptides, input_alist, input_mhc_seq,
         data_label, model, model_tag, results_dir,
         save_all_pred=False, input_pep_hla=False,
         encode9AA=False, neoantigen=False):

    start = datetime.now()

    config = str(data_label) + '_' + str(model_tag)
    model = str(model)

    # -----------------------------------------------------------------
    # read in data
    # -----------------------------------------------------------------

    print('read in data')
    test_predX = pd.read_csv(input_peptides, sep='\s+', header=0, engine='python')

    print('data dimension:')
    print(test_predX.shape)

    print('test_predX[0:2,]:')
    print(test_predX.iloc[0:2, ])

    # -----------------------------------------------------------------
    # read in HLA allele information for each cell line
    # -----------------------------------------------------------------

    allele_list = pd.read_csv(input_alist, sep='\s+', header=None, engine='python')
    print('HLA allele information:')
    print(allele_list.shape)
    print(allele_list.iloc[0:3, ])
    print(allele_list.iloc[1, 1].split(","))

    allele_dic = {}

    for i in allele_list.index:
        allele_dic[allele_list.iloc[i, 0]] = allele_list.iloc[i, 1].split(",")

    first2pairs = {k: allele_dic[k] for k in list(allele_dic)[:2]}
    print('first two entries of allele dictionary:')
    print(first2pairs)

    # -----------------------------------------------------------------
    # read in HLA sequences
    # -----------------------------------------------------------------

    hla_seq = pd.read_csv(input_mhc_seq, sep='\s+', header=None, engine='python')
    print('hla sequence information:')
    print(hla_seq.shape)
    print(hla_seq.iloc[0:3, :])

    # -----------------------------------------------------------------
    # encode all the HLA seqeunces
    # -----------------------------------------------------------------

    hla_encoding = {}

    for i in hla_seq.index:
        hla_name = hla_seq.iloc[i, 0]
        hla_aa = hla_seq.iloc[i, 1]
        hla_encoding[hla_name] = encode_peptide(hla_aa)

    print('finished hla_encode')

    # -----------------------------------------------------------------
    # For final evaluation of model: using input_test_pred
    # encode all possible HLA alleles for either binder or non-binder
    # -----------------------------------------------------------------

    test_idx = list()
    test_hla = list()
    test_encode = list()
    test_encode_hla = list()

    if input_pep_hla:
        print('prepare pep and hla as separate input')

    for i in test_predX.index:
        pep_encode = encode_peptide(test_predX.iloc[i, 0])
        sample_id = test_predX.iloc[i, 2]
        hla_allele = allele_dic[sample_id]

        for hla1 in hla_allele:
            if hla1 in hla_encoding:
                hla_encode = hla_encoding[hla1]
                test_idx.append(i)
                test_hla.append(hla1)
                if input_pep_hla:
                    test_encode.append(pep_encode.flatten())
                    test_encode_hla.append(hla_encode.flatten())
                else:
                    test_encode.append(np.concatenate((pep_encode, hla_encode)).flatten())
            else:
                print(f"{hla1} is not a key in the HLA dictionary.")

    df_test = test_predX.loc[test_idx, :]

    if neoantigen:
        df_test.columns.values[:3] = np.array(['peptide', 'key', 'sample'], dtype=object)
    else:
        df_test.columns.values[:3] = np.array(['peptide', 'y_true', 'sample'], dtype=object)

    test_encode = np.asarray(test_encode)
    test_encode_hla = np.asarray(test_encode_hla)

    end_encode = datetime.now()

    # -----------------------------------------------------------------
    # load model for prediction
    # -----------------------------------------------------------------

    model1 = tf.keras.models.load_model(model)

    # -----------------------------------------------------------------
    # Predict for weighted and unweighted models
    # -----------------------------------------------------------------
    if input_pep_hla:
        test_pred = model1.predict([test_encode, test_encode_hla])
    else:
        test_pred = model1.predict(test_encode)

    df_test['y_pred'] = test_pred

    # -----------------------------------------------------------------
    # Take maximum for peptide/sample and calculate classification
    # -----------------------------------------------------------------

    print("first 10 rows of prediction results")
    print(df_test.iloc[0:10, ])

    if encode9AA:
        var2grp = ['peptide_orig', 'y_true', 'sample']
    elif neoantigen:
        var2grp = ['key', 'sample']
    else:
        var2grp = ['peptide', 'y_true', 'sample']

    df_max = df_test.groupby(var2grp, as_index=False)['y_pred'].max()

    if not neoantigen:
        auc_roc = roc_auc_score(df_max["y_true"], df_max["y_pred"])
        pr_roc = average_precision_score(df_max["y_true"], df_max["y_pred"])

        print('evaluation (of max if MA)')
        print('auc ROC, auc PR')
        print([round(auc_roc, 3), round(pr_roc, 3)])

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

    fnm_max = config + '_test_pred.txt'
    df_max.to_csv(os.path.join(results_dir, fnm_max),
                  sep='\t', index=False, float_format='%.3e')
    print('Done with predictions!')

    if save_all_pred:
        fnm = config + '_test_pred_all.txt'
        df_test.to_csv(os.path.join(results_dir, fnm),
                       sep='\t', index=False, float_format='%.3e')
        subprocess.run(["gzip", os.path.join(results_dir, fnm)])

    end_pred = datetime.now()
    print("Start training =", start)
    print("Done encoding =", end_encode)
    print("Done =", end_pred)

    print('Finished all predictions. Bye!')


# -----------------------------------------------------------------
# parameters
# -----------------------------------------------------------------

parser = argparse.ArgumentParser(
    description='a pan-specific method for antigen presentation using deep learning.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "--input_peptides",
    type=str,
    dest="input_peptides",
    help="input: file name of peptides. Delimited by space or tab"
)
parser.add_argument(
    "--input_alist",
    type=str,
    dest="input_alist",
    help="input: file name for list of HLA-I alleles for all the samples. Delimited by space or tab"
)
parser.add_argument(
    "--input_mhc_seq",
    type=str,
    dest="input_mhc_seq",
    default="../data/NetMHCpan4_1_train/MHC_pseudo.dat",
    help="input: file of pseudo-sequence of HLA alleles. Delimited by space or tab"
)
parser.add_argument(
    "--data_label",
    type=str,
    dest="data_label",
    help="data label, which will be part of the output file names."
)
parser.add_argument(
    "--model",
    type=str,
    dest="model",
    help="input: file name of an .h5 file that saves a NN model."
)
parser.add_argument(
    "--model_tag",
    type=str,
    dest="model_tag",
    help="short string identifying the model."
)
parser.add_argument(
    "--results_dir", "-R",
    type=str,
    default="../results/test_data",
    help="path to directory where results are saved"
)
parser.add_argument(
    "--save_all_pred",
    dest="save_all_pred",
    action="store_true",
    default=False,
    help="option to save all predictions for all combinations of peptide and HLA in cell line"
)
parser.add_argument(
    "--input_pep_hla",
    dest="input_pep_hla",
    action="store_true",
    default=False,
    help="option to use an older version of the model where peptide and hla are separate inputs."
)
parser.add_argument(
    "--encode9AA",
    dest="encode9AA",
    action="store_true",
    default=False,
    help="whether the peptides are encoded as 9AA sequence otherwise they should be 15AA sequence."
)
parser.add_argument(
    "--neoantigen",
    dest="neoantigen",
    action="store_true",
    default=False,
    help="whether the testing data are for neoantigen prediction."
)

if __name__ == '__main__':
    arguments = parser.parse_args()
    main(**vars(arguments))
