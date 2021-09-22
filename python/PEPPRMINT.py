
from __future__ import absolute_import, division, print_function, unicode_literals
from __future__ import unicode_literals

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
# utility function to plot losses
# -------------------------------------------------------------------------

def plot_loss(m1, start, plot_path):
    # plot the training loss and accuracy
    end = len(m1.history['loss'])
    N   = np.arange(start, end)
    s   = slice(start,end)
    
    plt.style.use("ggplot")
    plt.figure(figsize=(4, 3), dpi=300)
    
    plt.plot(N, (m1.history["loss"][s]), label="train_loss")
    plt.plot(N, (m1.history["val_loss"][s]), label="val_loss")
    
    plt.xlabel("Epoch #")
    plt.ylabel("Loss")
    plt.legend()
    plt.subplots_adjust(left=0.2, right=0.98, top=0.98, bottom=0.2)
    plt.savefig(plot_path)
    plt.close()

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

# input_train  = "train_v4_el_multi_HLA_0.txt.gz"
# input_valid  = "train_v4_el_multi_HLA_1.txt.gz"
# input_test   = "../pMHCpan_data/test_v4_el_single_HLA.txt.gz"

# input_alist  = "allelelist.txt"
# input_MHCseq = "MHC_pseudo.txt"

# init_model   = "pMHCpan_800_bs_32_lr_0.001_e_50_layer_1_split0_one_input_Jun29_wsun_best_valid.h5"
# olabel       = "pMHC_balanced_0728_wsun"
# data_dir     = "../data/mixPep_data"
# info_dir     = "../data/NetMHCpan_train"
# fig_dir      = "../figures/mixPep"
# results_dir  = "../results/mixPep"

# n_layers     = 1
# hidden_size1 = 800
# hidden_size2 = 400
# n_epochs     = 2
# batch_size   = 32
# learn_rate   = 0.001
# n_iter       = 10
# decr_iter    = 2

# binder_weight = 1
# converge_e    = 1e-5
# converge_iter = 5
# use_class_weight    = False
# new_pi_weight       = 0.5
# save_all_iterations = True
# save_all_pred       = True
# save_model          = True

def main(input_train, input_valid, input_test, input_alist,
    input_MHCseq, init_model, olabel, 
    data_dir, info_dir, fig_dir, results_dir,
    n_layers, hidden_size1, hidden_size2,
    n_epochs, batch_size, learn_rate, n_iter, decr_iter, 
    converge_e = 1e-5, converge_iter = 5, 
    binder_weight = 10, new_pi_weight = 0.5, drop_rate = 0.5, 
    use_class_weight = False, 
    save_all_iterations = False, save_all_pred = False, 
    save_model = False):

    saved_args = locals()
    print("saved_args is", saved_args)

    # minimum prediction of presentation. This is needed when we calculate
    # weighted average to avoid divided by 0.
    start = datetime.now()

    MIN_PRED = 1e-7
    olabel   = str(olabel)

    # convergence checking
    converge_n  = 0
    decr_like_n = 0

    # -----------------------------------------------------------------
    # set model configuration
    # -----------------------------------------------------------------

    if n_layers > 2 :
        print('n_layers can be only 1 or 2:')
        raise ValueError

    hidden_sizes = [hidden_size1, hidden_size2]

    config = 'MA_' + str(hidden_sizes[0])

    if n_layers == 2:
        config = config + '_' + str(hidden_sizes[1])

    config = config + '_bs_' + str(batch_size) + '_lr_' + str(learn_rate)
    config = config + '_e_' + str(n_epochs) + '_layer_' + str(n_layers) + '_dropout_' + str(drop_rate)
    #config = config + 'decr' + str(decr_iter)
    config = config + '_new_pi_weight_' + str(new_pi_weight) + '_decriter_' + str(decr_iter)

    if use_class_weight:
        config = config + 'used_class_weight'

    config = config + '_' + olabel

    config

    # -----------------------------------------------------------------
    # read in data
    # Note: data should be pre-processed to not have duplicated peptides
    #   (i.e. same peptide binding to multiple cell lines)
    # -----------------------------------------------------------------

    file_train = os.path.join(data_dir, input_train)
    trainX     = pd.read_csv(file_train, sep='\t', header=0)

    file_valid = os.path.join(data_dir, input_valid)
    validX     = pd.read_csv(file_valid, sep='\t', header=0)

    file_test  = os.path.join(data_dir, input_test)
    testX      = pd.read_csv(file_test, sep='\t', header=0)

    print('training, validation, and testing data dimension:')
    print(trainX.shape)
    print(validX.shape)
    print(testX.shape)

    print('trainX[0:2,]:')
    print(trainX.iloc[0:2,])

    print('validX[0:2,]:')
    print(validX.iloc[0:2,])

    print('testX[0:2,]:')
    print(testX.iloc[0:2,])

    # -----------------------------------------------------------------
    # read in HLA allele information for each cell line
    # -----------------------------------------------------------------

    file_allele = os.path.join(info_dir, input_alist)
    allele_list = pd.read_csv(file_allele, sep=' ', header=None)
    allele_list.shape
    allele_list.iloc[0:3,]

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
    hla_seq.iloc[0:3,:]

    # -----------------------------------------------------------------
    # encode all the HLA seqeunces
    # -----------------------------------------------------------------

    hla_encoding = {}

    for i in hla_seq.index:
        hla_name = hla_seq.iloc[i,0]
        hla_aa   = hla_seq.iloc[i,1]
        hla_encoding[hla_name] = encode_peptide(hla_aa)

    # -----------------------------------------------------------------
    # encode peptide and hla data for model fitting for training data
    # encode all possible HLA alleles for binder
    # randomly chose an HLA allele for a non-binder
    # -----------------------------------------------------------------

    train_idx     = list()
    train_hla     = list()
    train_encode  = list()

    random.seed(2021)

    for i in trainX.index:
        pep1       = trainX.iloc[i,0]
        binder     = trainX.iloc[i,1]
        cell_line  = trainX.iloc[i,2]
        pep_encode = encode_peptide(pep1)
        hla_allele = allele_dic[cell_line]
        
        if binder:
            for hla1 in hla_allele:
                hla_encode = hla_encoding[hla1]
                train_idx.append(i)
                train_hla.append(hla1)
                train_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())
        else:
            hla1 = random.choice(hla_allele)
            hla_encode = hla_encoding[hla1]
            train_idx.append(i)
            train_hla.append(hla1)
            train_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())

    df_train = trainX.iloc[train_idx, [2,0,1]]
    df_train.rename(columns={"cell_line": "sample", "binder": "y_true"}, 
        inplace=True)
    df_train.insert(loc=1, column='hla', value=train_hla)
    df_train.reset_index(drop=True, inplace=True)

    print('df_train= ')
    print(df_train.iloc[0:6,])

    train_encode  = np.array(train_encode)
    train_encode.shape
    print('train_encode[0:2,]:')
    print(train_encode[0:6,])

    # -----------------------------------------------------------------
    # encode peptide and hla data for model fitting for validation
    # data : encode all possible HLA alleles for binder
    #        randomly chose an HLA allele for a non-binder
    # -----------------------------------------------------------------

    valid_idx     = list()
    valid_hla     = list()
    valid_encode  = list()

    random.seed(2021)

    for i in validX.index:
        pep1       = validX.iloc[i,0]
        binder     = validX.iloc[i,1]
        cell_line  = validX.iloc[i,2]
        pep_encode = encode_peptide(pep1)
        hla_allele = allele_dic[cell_line]
        
        if binder:
            for hla1 in hla_allele:
                hla_encode = hla_encoding[hla1]
                valid_idx.append(i)
                valid_hla.append(hla1)
                valid_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())
        else:
            hla1 = random.choice(hla_allele)
            hla_encode = hla_encoding[hla1]
            valid_idx.append(i)
            valid_hla.append(hla1)
            valid_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())

    df_valid  = validX.iloc[valid_idx,[2,0,1]]
    df_valid.rename(columns={"cell_line": "sample", "binder": "y_true"}, inplace=True)
    df_valid.insert(loc=1, column='hla', value=valid_hla)
    df_valid.reset_index(drop=True, inplace=True)

    print('df_valid= ')
    print(df_valid.iloc[0:4,])

    valid_encode  = np.array(valid_encode)
    valid_encode.shape

    print('valid_encode[0:2,]:')
    print(valid_encode[0:4,])

    # -----------------------------------------------------------------
    # encode peptide and hla data for for testing
    #   encode all possible HLA alleles for either binder or non-binder
    # -----------------------------------------------------------------

    test_idx     = list()
    test_hla     = list()
    test_encode  = list()

    for i in testX.index:
        pep_encode = encode_peptide(testX.iloc[i,0])
        binder     = testX.iloc[i,1]
        cell_line  = testX.iloc[i,2]
        hla_allele = allele_dic[cell_line]
        
        for hla1 in hla_allele:
            hla_encode = hla_encoding[hla1]
            test_idx.append(i)
            test_hla.append(hla1)
            test_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())

    df_test  = testX.iloc[test_idx,[2,0,1]]
    # ad-hoc here to change the header of testing data from HLA to sample
    df_test.rename(columns={"HLA": "sample", "binder": "y_true"}, inplace=True)
    #df_test.rename(columns={"cell_line": "sample", "binder": "y_true"}, inplace=True)
    df_test.insert(loc=1, column='hla', value=test_hla)
    df_test.reset_index(drop=True, inplace=True)

    print('df_test= ')
    print(df_test.iloc[0:2,])

    test_encode  = np.array(test_encode)
    test_encode.shape
    print('test_encode[0:2,]:')
    print(test_encode[0:2,])

    end_encode = datetime.now()

    # -----------------------------------------------------------------
    # load initial model trained using single allele data
    # -----------------------------------------------------------------

    file_h5 = os.path.join(results_dir, init_model)
    pMHCpan = tf.keras.models.load_model(file_h5)

    # -----------------------------------------------------------------
    # check the predicted presentation probability using initial model
    # -----------------------------------------------------------------

    y_pred_train = pMHCpan.predict(train_encode)
    y_pred_valid = pMHCpan.predict(valid_encode)

    df_train['y_pred_pMHC']  = y_pred_train
    df_valid['y_pred_pMHC']  = y_pred_valid
    df_test['y_pred_pMHC']   = pMHCpan.predict(test_encode)

    y_pred_train.max()
    y_pred_train.min()

    y_pred_train[y_pred_train < MIN_PRED]   = MIN_PRED
    y_pred_train[y_pred_train > 1-MIN_PRED] = 1-MIN_PRED

    y_pred_valid[y_pred_valid < MIN_PRED]   = MIN_PRED
    y_pred_valid[y_pred_valid > 1-MIN_PRED] = 1-MIN_PRED

    df_train['y_pred']  = y_pred_train
    df_valid['y_pred']  = y_pred_valid

    max1 = df_train.groupby(['sample','peptide','y_true'], as_index=False)['y_pred'].max()
    max1.groupby('y_true', as_index=False)['y_pred'].describe()

    max1 = df_valid.groupby(['sample','peptide','y_true'], as_index=False)['y_pred'].max()
    max1.groupby('y_true', as_index=False)['y_pred'].describe()

    # -----------------------------------------------------------------
    # normalize the prediction within each group so that
    # the summation of weights is 1
    # -----------------------------------------------------------------

    y1_train = (df_train['y_true']==1)
    y1_valid = (df_valid['y_true']==1)

    y0_train = (df_train['y_true']==0)
    y0_valid = (df_valid['y_true']==0)

    wt_train = (df_train.loc[y1_train,:]).groupby(['sample','peptide'])['y_pred'].transform(lambda x: x / x.sum())
    wt_valid = (df_valid.loc[y1_valid,:]).groupby(['sample','peptide'])['y_pred'].transform(lambda x: x / x.sum())

    df_train.loc[y1_train, 'p_ijk'] = wt_train
    df_valid.loc[y1_valid, 'p_ijk'] = wt_valid

    df_train.loc[y0_train, 'p_ijk'] = 1
    df_valid.loc[y0_valid, 'p_ijk'] = 1

    print('initilized weights')
    print(df_train.iloc[0:6,])

    # double check the summation of weights are indeed 1
    sum1 = df_train.groupby(['sample','peptide','y_true'], as_index=False)['p_ijk'].sum()
    print('check that the weights sum up to 1')
    print(sum1.groupby('y_true')['p_ijk'].describe())

    sum1 = df_valid.groupby(['sample','peptide','y_true'], as_index=False)['p_ijk'].sum()
    sum1.groupby('y_true')['p_ijk'].describe()

    # -----------------------------------------------------------------
    # set up the NN model
    # -----------------------------------------------------------------
    #define class weights
    bind_class_weight = {
        0: 1.,
        1: binder_weight}

    if use_class_weight:
        print("binder class weight used in NN ")
        print(bind_class_weight)

    # NN
    input1 = keras.Input(shape=(train_encode.shape[1],), name='input')
    hidden = layers.Dense(hidden_sizes[0], activation='relu', name='hidden0')(input1)
    hidden = layers.Dropout(rate=drop_rate, name='dropout_hidden0')(hidden)

    if n_layers > 1:
        hidden = layers.Dense(hidden_sizes[1], activation='relu', name='hidden1')(hidden)
        hidden = layers.Dropout(rate=drop_rate, name='dropout_hidden1')(hidden)

    if n_layers > 2:
        hidden = layers.Dense(hidden_sizes[2], activation='relu', name='hidden2')(hidden)
        hidden = layers.Dropout(rate=drop_rate, name='dropout_hidden2')(hidden)

    prob = layers.Dense(1, activation='sigmoid', name='output')(hidden)

    mixPep = keras.Model(input1, prob)
    mixPep.summary()

    adam1 = keras.optimizers.Adam(lr=learn_rate)
    mixPep.compile(optimizer=adam1, loss='binary_crossentropy', metrics=['accuracy'])

    # =================================================================
    # start the EM algorithm
    # =================================================================
    start_EM = datetime.now()

    logLiks = np.zeros(n_iter)

    for idx in range(n_iter):
        print('----------------------------------------------------------')
        print('iteration ' + str(idx+1))
        print('----------------------------------------------------------')

        # -----------------------------------------------------------------
        # The M step 1: re-estimate \pi_{ik}
        # -----------------------------------------------------------------

        y1_train = (df_train['y_true']==1)
        grp_pi = (df_train.loc[y1_train,:]).groupby(['sample','hla'], as_index=False)
        df_pi  = grp_pi['p_ijk'].sum()
        # print(df_pi)

        # as a common issue for mixture model, some component may be empty.
        # this will create numerical issues and here we set a lower-bound.
        df_pi.loc[df_pi['p_ijk'] < MIN_PRED, 'p_ijk'] = MIN_PRED

        grp_sample = df_pi.groupby('sample', as_index=False)
        pi = grp_sample['p_ijk'].transform(lambda x: x / x.sum())
        df_pi.loc[:,'pi'] = np.array(pi)
        df_pi

        if idx == 0:
            df_pi0 = df_pi

        df_pi['pi'] = new_pi_weight*df_pi['pi'] + (1 - new_pi_weight)*df_pi0['pi']

        # -----------------------------------------------------------------
        # the M step 2: estimate peptide presentation model
        # -----------------------------------------------------------------

        train_y_true = np.array(df_train['y_true'])
        valid_y_true = np.array(df_valid['y_true'])

        wt_train = np.array(df_train['p_ijk'])
        wt_valid = np.array(df_valid['p_ijk'])

        if use_class_weight:
            print('fit NN with class weights')
            m0 = mixPep.fit(train_encode, train_y_true, epochs=n_epochs,
                batch_size=batch_size, verbose=2, sample_weight=wt_train,
                class_weight=bind_class_weight,
                validation_data=(valid_encode, valid_y_true, wt_valid))
        else:
            print('fit WITHOUT class weights')
            m0 = mixPep.fit(train_encode, train_y_true, epochs=n_epochs,
                batch_size=batch_size, verbose=2, sample_weight=wt_train,
                validation_data=(valid_encode, valid_y_true, wt_valid))

        # -----------------------------------------------------------------
        # The E step, estimate probablities across mixtures
        # -----------------------------------------------------------------

        y_pred_train = mixPep.predict(train_encode)
        y_pred_valid = mixPep.predict(valid_encode)

        # print('predicted train (E-step)')
        # print(y_pred_train[0:6])
        #print(max(y_pred_train))

        y_pred_train[y_pred_train < MIN_PRED]   = MIN_PRED
        y_pred_train[y_pred_train > 1-MIN_PRED] = 1-MIN_PRED

        y_pred_valid[y_pred_valid < MIN_PRED]   = MIN_PRED
        y_pred_valid[y_pred_valid > 1-MIN_PRED] = 1-MIN_PRED

        df_train['y_pred'] = y_pred_train
        df_valid['y_pred'] = y_pred_valid
        #print(max(df_train['y_pred']))
        #print(min(df_train['y_pred']))
        #print(df_train.groupby(['y_true', 'hla'], as_index=False)['y_pred'].summary())

        df_train['wt_y_pred'] = y_pred_train
        df_valid['wt_y_pred'] = y_pred_valid

        # -----------------------------------------------------------------
        # multiply predicted density by estimates of \pi
        # -----------------------------------------------------------------

        for ix in df_pi.index:
            #print(df_pi)
            #print(df_pi0)
            sample1 = df_pi.loc[ix,'sample']
            hla1    = df_pi.loc[ix,'hla']
            pi1     = df_pi.loc[ix,'pi']
            
            wUpdate = (df_train['sample']==sample1) & (df_train['hla']==hla1) & (df_train['y_true']==1)
            df_train.loc[wUpdate,'wt_y_pred'] = df_train.loc[wUpdate,'wt_y_pred']*pi1
            
            wUpdate = (df_valid['sample']==sample1) & (df_valid['hla']==hla1) & (df_valid['y_true']==1)
            df_valid.loc[wUpdate,'wt_y_pred'] = df_valid.loc[wUpdate,'wt_y_pred']*pi1

        y1_train = (df_train['y_true']==1)
        y1_valid = (df_valid['y_true']==1)

        # normalize the summation to one for each sample and peptide, across HLA alleles
        wt_train = (df_train.loc[y1_train,:]).groupby(['sample','peptide'])['wt_y_pred'].transform(lambda x: x / x.sum())
        wt_valid = (df_valid.loc[y1_valid,:]).groupby(['sample','peptide'])['wt_y_pred'].transform(lambda x: x / x.sum())

        df_train.loc[y1_train, 'p_ijk'] = wt_train
        df_valid.loc[y1_valid, 'p_ijk'] = wt_valid

        # print(df_train.groupby(['y_true'], as_index=False)['p_ijk'].describe())
        # print(df_valid.groupby(['y_true'], as_index=False)['p_ijk'].describe())


        # -----------------------------------------------------------------
        # calculate the log-likelihood
        # -----------------------------------------------------------------

        df_train_sum = df_train.groupby(['sample','peptide','y_true'], as_index=False)['wt_y_pred'].sum()
        df_valid_sum = df_valid.groupby(['sample','peptide','y_true'], as_index=False)['wt_y_pred'].sum()

        df_train_sum.groupby(['y_true'], as_index=False)['wt_y_pred'].describe()

        print('df_train_sum')
        print(df_train_sum)


        y0_train = (df_train_sum['y_true']==0)
        y0_valid = (df_valid_sum['y_true']==0)

        y1_train = (df_train_sum['y_true']==1)
        y1_valid = (df_valid_sum['y_true']==1)

        if use_class_weight:
            log_lik_train = (np.log(binder_weight* df_train_sum.loc[y1_train,'wt_y_pred'])).sum()
            log_lik_valid = (np.log(binder_weight* df_valid_sum.loc[y1_valid,'wt_y_pred'])).sum()
        else:
            log_lik_train = (np.log(df_train_sum.loc[y1_train,'wt_y_pred'])).sum()
            log_lik_valid = (np.log(df_valid_sum.loc[y1_valid,'wt_y_pred'])).sum()

        log_lik_train = log_lik_train + (np.log(1 - df_train_sum.loc[y0_train,'wt_y_pred'])).sum()
        log_lik_valid = log_lik_valid + (np.log(1 - df_valid_sum.loc[y0_valid,'wt_y_pred'])).sum()

        print('')
        print('~~~~~~ log_lik_train and log_lik_valid ~~~~~~')
        print('')
        print(str(log_lik_train) + '   ' + str(log_lik_valid))
        print('')

        logLiks[idx] = log_lik_train

        # if likelihood decrease, exit the loop.
        if idx > 0:
            if(log_lik_train + 1e-5*np.abs(log_lik_train) < logLiks[idx-1]):
                warnings.warn('train log likelihood decreased.')
                print('train log decreased in iteration =')
                print(idx+1)
                decr_like_n = decr_like_n+1
                
                print('~~~~ total number of log likelihood decreases = ')
                print(decr_like_n)
                
                if(decr_like_n > decr_iter):
                    print('**** EARLY STOPPING: log likelihood decreased decr_iter times')
                    break
            else: 
                decr_like_n = 0

        # save the df_w and df_pi for output at the end of the loop
        # in case the likelihood decreases and we need to remember them
        df_pi0 = df_pi

        # -----------------------------------------------------------------
        # Evaluate results using Internal encoding for TRAINING dataset:
        #     for each peptide binder, take the HLA with maximum prediction
        #     as the paired HLA, and evaluate results
        # -----------------------------------------------------------------

        df_train_grp    = df_train.groupby(['sample','peptide','y_true'], as_index=False)
        df_train_pMHC   = df_train_grp['y_pred_pMHC'].max()
        df_train_mixPep = df_train_grp['y_pred'].max()

        auc_train = roc_auc_score(df_train_pMHC["y_true"], df_train_pMHC["y_pred_pMHC"])
        pr_train  = average_precision_score(df_train_pMHC["y_true"], df_train_pMHC["y_pred_pMHC"])

        auc_train_mix = roc_auc_score(df_train_mixPep["y_true"], df_train_mixPep["y_pred"])
        pr_train_mix  = average_precision_score(df_train_mixPep["y_true"], df_train_mixPep["y_pred"])

        print("----------------------------------------------------")
        print('Max pred with INTERNAL encoding of TRAIN SET')
        print("----------------------------------------------------")
        print('AUC ROC [pMHC, mixPep]')
        print([round(auc_train,3), round(auc_train_mix,3)])

        print('AUC PR [pMHC,mixPep]')
        print([round(pr_train,3), round(pr_train_mix, 3)])

        # -----------------------------------------------------------------
        # Evaluate results using INTERNAL encoding for VALIDATION set:
        #     for each peptide binder, take the HLA with maximum prediction
        #     as the paired HLA, and evaluate results
        # -----------------------------------------------------------------

        df_valid_grp    = df_valid.groupby(['sample', 'peptide','y_true'], as_index=False)
        df_valid_pMHC   = df_valid_grp['y_pred_pMHC'].max()
        df_valid_mixPep = df_valid_grp['y_pred'].max()

        auc_valid = roc_auc_score(df_valid_pMHC["y_true"], df_valid_pMHC["y_pred_pMHC"])
        pr_valid  = average_precision_score(df_valid_pMHC["y_true"], df_valid_pMHC["y_pred_pMHC"])

        auc_valid_mix = roc_auc_score(df_valid_mixPep["y_true"], df_valid_mixPep["y_pred"])
        pr_valid_mix  = average_precision_score(df_valid_mixPep["y_true"], df_valid_mixPep["y_pred"])

        print("----------------------------------------------------")
        print('Max pred with Internal encoding of VALIDATION SET')
        print("----------------------------------------------------")
        print('AUC ROC [pMHC, mixPep]')
        print([round(auc_valid,3), round(auc_valid_mix,3)])

        print('AUC PR [pMHC,mixPep]')
        print([round(pr_valid,3), round(pr_valid_mix, 3)])

        # -----------------------------------------------------------------
        # Evaluate results using External for TEST set:
        #     for each peptide binder, take the HLA with maximum prediction
        #     as the paired HLA, and evaluate results
        # -----------------------------------------------------------------

        df_test['y_pred']  = mixPep.predict(test_encode)
        df_test

        df_test_grp    = df_test.groupby(['sample', 'peptide','y_true'], as_index=False)
        df_test_pMHC   = df_test_grp['y_pred_pMHC'].max()
        df_test_mixPep = df_test_grp['y_pred'].max()

        auc_test = roc_auc_score(df_test_pMHC["y_true"], df_test_pMHC["y_pred_pMHC"])
        pr_test  = average_precision_score(df_test_pMHC["y_true"], df_test_pMHC["y_pred_pMHC"])

        auc_test_mix = roc_auc_score(df_test_mixPep["y_true"], df_test_mixPep["y_pred"])
        pr_test_mix  = average_precision_score(df_test_mixPep["y_true"], df_test_mixPep["y_pred"])

        print("----------------------------------------------------")
        print('Max pred with external encoding of TEST SET')
        print("----------------------------------------------------")
        print('AUC ROC [pMHC, mixPep]')
        print([round(auc_test,3), round(auc_test_mix,3)])

        print('AUC PR [pMHC,mixPep]')
        print([round(pr_test,3), round(pr_test_mix, 3)])

        if save_all_iterations:
            
            if save_all_pred:
                fnm  = config + '_iter'+ str(idx+1) +'_test_pred_all.txt'
                df_test.to_csv(os.path.join(results_dir, fnm),
                    sep='\t', index=False, float_format='%.3e')
                subprocess.run(["gzip", os.path.join(results_dir, fnm)])
                
                fnm  = config + '_iter'+ str(idx+1) +'_valid_pred_all.txt'
                df_valid.to_csv(os.path.join(results_dir, fnm),
                    sep='\t', index=False, float_format='%.3e')
                subprocess.run(["gzip", os.path.join(results_dir, fnm)])
                
                fnm  = config + '_iter' + str(idx+1) +'_train_pred_all.txt'
                df_train.to_csv(os.path.join(results_dir, fnm),
                    sep='\t', index=False, float_format='%.3e')
                subprocess.run(["gzip", os.path.join(results_dir, fnm)])
                
                fnm_pi =  config +  'iter' + str(idx+1) +'_pi.txt'
                df_pi0.to_csv(os.path.join(results_dir, fnm_pi),
                    sep='\t', index=False, float_format='%.3e')
                
            if save_model:
                fnm = config + '_iter' + str(idx+1) + '.h5'
                model_path  = os.path.join(results_dir, fnm)
                mixPep.save(model_path)

        # -----------------------------------------------------------------
        # Checks for early stopping if the likelihood increases is small
        # for converge_itr consecutative iterations.
        # -----------------------------------------------------------------

        if idx > 0:
            if(log_lik_train < logLiks[idx-1] + converge_e*np.abs(logLiks[idx-1])):
                converge_n = converge_n + 1
            else:
                converge_n = 0
            
            if(converge_n >= converge_iter):
                print('**** EARLY STOPPING: convergence of loglikelihood ****')
                break
            
    end_EM = datetime.now()

    # -----------------------------------------------------------------
    # Save predictions: Assuming log likelihood did not decrease
    # -----------------------------------------------------------------

    if (not save_all_iterations) and save_all_pred:
        fnm  = config + '_iter'+ str(idx+1) +'_test_pred_all.txt'
        df_test.to_csv(os.path.join(results_dir, fnm),
            sep='\t', index=False, float_format='%.3e')
        subprocess.run(["gzip", os.path.join(results_dir, fnm)])
        
        fnm  = config + '_iter'+ str(idx+1) +'_valid_pred_all.txt'
        df_valid.to_csv(os.path.join(results_dir, fnm),
            sep='\t', index=False, float_format='%.3e')
        subprocess.run(["gzip", os.path.join(results_dir, fnm)])
        
        fnm  = config + '_iter' + str(idx+1) +'_train_pred_all.txt'
        df_train.to_csv(os.path.join(results_dir, fnm),
            sep='\t', index=False, float_format='%.3e')
        subprocess.run(["gzip", os.path.join(results_dir, fnm)])
        
        fnm_pi =  config +  '_iter' + str(idx+1) +'_pi.txt'
        df_pi0.to_csv(os.path.join(results_dir, fnm_pi),
            sep='\t', index=False, float_format='%.3e')
        
    if (not save_all_iterations) and save_model:
        fnm = config + '_iter' + str(idx+1) + '.h5'
        model_path  = os.path.join(results_dir, fnm)
        mixPep.save(model_path)

    print('log_lik_train across iterations:')
    print(logLiks)
    print("Start training =", start)
    print("Done encoding =", end_encode)
    print("Start EM =", start_EM)
    print("Done EM =", end_EM)
    print('Finished all iterations of training. Bye!')

# -----------------------------------------------------------------
# parameters
# -----------------------------------------------------------------

parser = argparse.ArgumentParser(
    description='a pan-specific method for antigen presentation using deep learning.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "--input_train",
    type = str,
    dest = "input_train",
    help = "input: file name of training data"
)
parser.add_argument(
    "--input_valid",
    type = str,
    dest = "input_valid",
    help = "input: file name of validation data"
)
parser.add_argument(
    "--input_test",
    type = str,
    dest = "input_test",
    help = "input: file name of testing data"
)
parser.add_argument(
    "--input_alist",
    type = str,
    dest = "input_alist",
    default = "allelelist.txt",
    help = "input: list of HLA alleles for all the samples"
)
parser.add_argument(
    "--input_MHCseq",
    type = str,
    dest = "input_MHCseq",
    default = "MHC_pseudo.txt",
    help = "input: sequence of different HLA alleles"
)
parser.add_argument(
    "--init_model",
    type = str,
    dest = "init_model",
    help = "input: model that initializes the weights"
)
parser.add_argument(
    "--data_dir", "-D",
    type = str,
    dest = "data_dir",
    default = "../data/MA_data",
    help = "directory where data are placed"
)
parser.add_argument(
    "--info_dir", "-I",
    type = str,
    dest = "info_dir",
    default = "../data/NetMHCpan4_1_train",
    help = "directory where HLA allele informaton are placed"
)
parser.add_argument(
    "--fig_dir", "-F",
    type = str,
    dest = "fig_dir",
    default = "../figures/PEPPRMINT",
    help = "directory where figures are saved"
)
parser.add_argument(
    "--results_dir", "-R",
    type = str,
    dest = "results_dir",
    default = "../results/PEPPRMINT",
    help = "directory where results are saved"
)
parser.add_argument(
    "--n_layers", "-L",
    type = int,
    dest = "n_layers",
    default = 1,
    help = "number of layers for encoder and decoder: 1 or 2"
)
parser.add_argument(
    "--hidden_size1",
    type = int,
    dest = "hidden_size1",
    default = 800,
    help = "number of nodes for hidden layer 1"
)
parser.add_argument(
    "--hidden_size2",
    type = int,
    dest = "hidden_size2",
    default = 400,
    help = "number of nodes for hidden layer 2"
)
parser.add_argument(
    "--batch_size", "-M",
    type = int,
    dest = "batch_size",
    default = 32,
    help = "batch size used when training"
)
parser.add_argument(
    "--learn_rate",
    type = float,
    dest = "learn_rate",
    default = 1e-3,
    help = "learning rate when training"
)
parser.add_argument(
    "--olabel",
    type = str,
    dest = "olabel",
    help = "Label to append to all output for identification"
)
parser.add_argument(
    "--n_iter",
    type = int,
    dest = "n_iter",
    default = 10,
    help = "number of iterations for which to train"
)
parser.add_argument(
    "--n_epochs", "-e",
    type = int,
    dest = "n_epochs",
    default = 5,
    help = "number of epochs for which to train per iteration"
)
parser.add_argument(
    "--decr_iter",
    type = int,
    dest = "decr_iter",
    default = 2,
    help = "number of time log like can decrease before stop training"
)
parser.add_argument(
    "--binder_weight",
    type = float,
    dest = "binder_weight",
    default = 10,
    help = "class weight for binders (vs. nonbinder). used only if flag 'use_class_weight'"
)
parser.add_argument(
    "--new_pi_weight",
    type = float,
    default = 0.5,
    dest = "new_pi_weight",
    help = "update pi by a weighted average of new estimate of pi \
    and the pi estimate from previous iteration. This is the weight \
    for the new estimate of pi"
)
parser.add_argument(
    "--drop_rate",
    type = float,
    default = 0.5,
    dest = "drop_rate",
    help = "Dropout rate used in NN layers"
)
parser.add_argument(
    "--use_class_weight",
    dest = "use_class_weight",
    action = "store_true",
    default = False,
    help = "Use binder class weights when fitting NN"
)
parser.add_argument(
    "--save_all_iterations",
    dest = "save_all_iterations",
    action = "store_true",
    default = False,
    help = "whether seek to save results for all iterations"
)
parser.add_argument(
    "--save_all_pred",
    dest = "save_all_pred",
    action = "store_true",
    default = False,
    help = "save all validation predictions"
)
parser.add_argument(
    "--save_model",
    dest = "save_model",
    action = "store_true",
    default = False,
    help = "save the model to a h5 file"
)



if __name__ == '__main__':
    arguments = parser.parse_args()
    main(**vars(arguments))

