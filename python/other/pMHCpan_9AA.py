
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

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Model

from tensorflow.keras.layers import Input, Dense, Flatten, Reshape
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

# input_train    = "train_v4_el_single_HLA_9AA_0.txt.gz"
# input_validate = "train_v4_el_single_HLA_9AA_1.txt.gz"
# olabel         = "split0_Aug5_9AA_wsun"
# data_dir       = "../data/pMHCpan_data"
# info_dir       = "../data/NetMHCpan_train"
# fig_dir        = "../figures/pMHCpan"
# results_dir    = "../results/pMHCpan"

# validating_frac = 0.2

# hidden_size1 = 800
# hidden_size2 = 400
# n_layers     = 1

# n_epochs     = 3
# patience     = 5

# batch_size    = 32
# learn_rate    = 1e-3
# binder_weight = 5
# n_iter        = 3

# use_class_weight = False
# use_nn_encoding  = False
# save_validate_pred  = True

def main(input_train, input_validate, olabel,
    data_dir, info_dir, fig_dir, results_dir, 
    validating_frac, n_layers,  hidden_size1, hidden_size2,  
    n_epochs, n_iter, batch_size, learn_rate, binder_weight = 5, 
    converge_e = 1e-5, converge_itr = 5, patience = 5, 
    use_class_weight = False, use_nn_encoding = False, 
    save_validate_pred = True):

    saved_args = locals()
    print("saved_args is", saved_args)

    start_time = datetime.now()
    print("start time :", start_time)

    MIN_PRED = 1e-7

    if not use_class_weight:
        binder_weight = 1.0

    # -----------------------------------------------------------------
    # set model configuration
    # -----------------------------------------------------------------
    olabel = str(olabel)

    if n_layers > 2 :
            print('n_layers can be only 1 or 2:')
            raise ValueError

    hidden_sizes = [hidden_size1, hidden_size2]

    config = 'pMHCpan_' + str(hidden_sizes[0])

    if n_layers == 2:
        config = config + '_' + str(hidden_sizes[1])

    config = config + '_bs_' + str(batch_size) + '_lr_' + str(learn_rate)
    config = config + '_e_' + str(n_epochs) + '_layer_' + str(n_layers) 

    if use_class_weight: 
        config = config + '_use_class_weight' + str(binder_weight)

    if use_nn_encoding:  
        config = config + '_use_nn_encoding'

    config = config + '_' + olabel

    print(config)
    print('')

    # specify the name of the model with best performance on validation data set
    # in the training process, for saving the model parameters in later code
    checkpoint_path = \
      os.path.join(results_dir, config + '_best_valid.h5')
      
    # specify the name of the file to save the last model when the training stops
    save_path = \
      os.path.join(results_dir, config + '_last_saved.h5')

    # -----------------------------------------------------------------
    # read in data
    # -----------------------------------------------------------------

    file_name = os.path.join(data_dir, input_train)
        
    print('input file:')
    print(file_name)
    print('')

    df = pd.read_csv(file_name, sep='\t', header=0)
    df.shape
    df.iloc[0:3,]

    # -----------------------------------------------------------------
    # split_training_validation_data
    # -----------------------------------------------------------------

    if input_validate:
        trainX     = df
        file_valid = os.path.join(data_dir, input_validate)
        validX     = pd.read_csv(file_valid, sep='\t', header=0)
        print('file name for validation data is provided, will not split data\n')
    else:
        print('file name for validation data is not provided, will split data\n')
        dfT = df.T
        trainX, validX = train_test_split(dfT, 
            test_size=validating_frac, random_state=1999)

    print('training and validation data dimension:')
    print(trainX.shape)
    print(validX.shape)

    print('trainX[0:2,]:')
    print(trainX.iloc[0:2,])
    print('validX[0:2,]:')
    print(validX.iloc[0:2,])


    # randomly shuffle training data
    random.seed(10)
    trainX = trainX.sample(frac=1.0)
    trainX.reset_index(drop=True, inplace=True)
    print('training data dimension:')
    print(trainX.shape)
    print('trainX[0:2,]:')
    print(trainX.iloc[0:2,])

    # randomly sample 30% of validation data
    random.seed(10)
    validX = validX.sample(frac=0.3)
    validX.reset_index(drop=True, inplace=True)
    print('validation data dimension:')
    print(validX.shape)
    print('validX[0:2,]:')
    print(validX.iloc[0:2,])

    # -----------------------------------------------------------------
    # read in HLA seqeunces
    # this input file is a modified version of the file MHC_pseudo.dat
    # from NNAlign_MA_testsuite, by replacing \t with ' ', and more than 
    # one space by single space. 
    # -----------------------------------------------------------------

    file_hla = os.path.join(info_dir, "MHC_pseudo.txt")
    hla_seq  = pd.read_csv(file_hla, sep=' ', header=None)
    hla_seq.shape
    hla_seq.iloc[0:3,0:5]

    # -----------------------------------------------------------------
    # encode all the HLA seqeunces
    # -----------------------------------------------------------------

    hla_encoding = {}

    for i in hla_seq.index:
        hla_name = hla_seq.iloc[i,0]
        hla_aa   = hla_seq.iloc[i,1]
        hla_encoding[hla_name] = encode_peptide(hla_aa)

    print('done HLA encoding')

    # -----------------------------------------------------------------
    # encoding peptide and hla data
    # -----------------------------------------------------------------

    train_encode = list()
    valid_encode = list()

    hla_train_encode = list()
    hla_valid_encode = list()

    for i in trainX.index:
        pep_encode = encode_peptide(trainX.loc[i,'pep_core'])
        hla_encode = hla_encoding[trainX.loc[i,'HLA']]
        if use_nn_encoding:
            train_encode.append((pep_encode).flatten())
            hla_train_encode.append((hla_encode).flatten())
        else:
            train_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())

    for i in validX.index:
        pep_encode = encode_peptide(validX.loc[i,'pep_core'])
        hla_encode = hla_encoding[validX.loc[i,'HLA']]
        if use_nn_encoding:
            valid_encode.append((pep_encode).flatten())
            hla_valid_encode.append((hla_encode).flatten())
        else:
            valid_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())

    train_encode = np.array(train_encode)
    valid_encode = np.array(valid_encode)

    hla_train_encode = np.array(hla_train_encode)
    hla_valid_encode = np.array(hla_valid_encode)


    print(train_encode.shape)
    print(hla_train_encode.shape)

    print(valid_encode.shape)
    print(hla_valid_encode.shape)

    print('train_encode[0:2,]:')
    print(train_encode[0:2,])

    if use_nn_encoding:
        print('hla_train_encode[0:2,]:')
        print(hla_train_encode[0:2,])

    print('valid_encode[0:2,]:')
    print(valid_encode[0:2,])
    print('')

    df_train  = trainX.iloc[:,[2,0,1,3,4,5,6]]
    df_train.rename(columns={"HLA": "hla", "binder": "y_true", "length": "pep_length"}, 
        inplace=True)

    print('df_train= ')
    print(df_train.iloc[0:2,])

    df_valid  = validX.iloc[:,[2,0,1,3,4,5,6]]
    df_valid.rename(columns={"HLA": "hla", "binder": "y_true", "length": "pep_length"}, 
        inplace=True)

    print('df_valid= ')
    print(df_valid.iloc[0:2,])

    # -----------------------------------------------------------------
    # model setup
    # -----------------------------------------------------------------

    if use_class_weight:
        print('')
        print("----------------------------------------------------")
        print("~~~~~ binder class weight used in NN ")
        print("----------------------------------------------------")
        print(binder_weight)
        print('')

    input1 = Input(shape=(train_encode.shape[1],), name='input_peptide')

    if use_nn_encoding:
        input2 = Input(shape=(hla_train_encode.shape[1],), name='input_hla')
        
        # Peptide only NN
        hidden1 = Dense(64, activation='relu', name='hidden1_pep')(input1)
        hidden1 = Dropout(rate=0.5, name='dropout_hidden1_pep')(hidden1)
        
        # HLA only NN
        hidden2 = Dense(64, activation='relu', name='hidden2_hla')(input2)
        hidden2 = Dropout(rate=0.5, name='dropout_hidden2_hla')(hidden2)
        
        # combine the two branches
        x = concatenate([hidden1, hidden2])
    else:
        x = input1

    # apply layers on the combined outputs
    x = Dense(hidden_sizes[0], activation='relu', name='all_hidden0')(x)
    x = Dropout(rate=0.5, name='all_dropout_hidden0')(x)

    prob = Dense(1, activation='sigmoid', name='output')(x)

    if use_nn_encoding:
        pMHCpan = Model(inputs=[input1, input2], outputs=prob)
    else:
        pMHCpan = Model(inputs=input1, outputs=prob)

    print(pMHCpan.summary())

    # -----------------------------------------------------------------
    # model configuration
    # -----------------------------------------------------------------

    adam1 = keras.optimizers.Adam(learning_rate=learn_rate)

    METRICS = [
        tf.keras.metrics.BinaryAccuracy(name = "accuracy"),
        tf.keras.metrics.AUC(name = "auc"), 
        keras.metrics.AUC(name='auprc', curve='PR')
    ]

    pMHCpan.compile(loss=tf.keras.losses.BinaryCrossentropy(), optimizer= adam1, 
        metrics=METRICS)

    # =================================================================
    # start the EM algorithm
    # =================================================================

    start_EM = datetime.now()

    logLiks = np.zeros(n_iter)

    for idx in range(n_iter):
        print('')
        print('----------------------------------------------------------')
        print('iteration ' + str(idx))
        print('----------------------------------------------------------')

        # -----------------------------------------------------------------
        # The M step 1: re-estimate \w_{kth}
        # -----------------------------------------------------------------

        # w: prob for starting position of binding core for a pep_length
        gby   = ['pep_length', 'start_pos']
        grp_w = df_train[df_train['y_true']==1].groupby(gby, as_index=False)
        df_w  = grp_w['weight'].sum()

        grp_hp = df_w.groupby('pep_length', as_index=False)
        w      = grp_hp['weight'].transform(lambda x: x / x.sum())
        df_w.loc[:,'w'] = np.array(w)

        print('df_w')
        print(df_w)

        # -----------------------------------------------------------------
        # the M step 2: estimate peptide presentation model
        #   evaluate the model using the validation set
        # -----------------------------------------------------------------

        train_y_true = np.array(df_train['y_true'])
        valid_y_true = np.array(df_valid['y_true'])

        w_train = np.array(df_train['weight'])
        w_valid = np.array(df_valid['weight'])

        if use_class_weight:
            print('')
            print("~~~~~ fit NN with class weights")
            class_weights = {0:1, 1:binder_weight}
        else: 
            print('')
            print("~~~~ fit NN without class weights")
            class_weights = {0:1, 1:1}

        callback = tf.keras.callbacks.EarlyStopping(monitor='val_accuracy', 
            patience=patience)

        # this check_point saves the model with best performance on validation data
        check_point = tf.keras.callbacks.ModelCheckpoint(
            checkpoint_path,
            monitor="val_accuracy",
            verbose=1,
            save_best_only=True,
            mode="auto"
        )

        if use_nn_encoding:
            m1 = pMHCpan.fit([train_encode, hla_train_encode], train_y_true, \
                epochs=n_epochs, batch_size=batch_size, verbose=2, sample_weight=w_train, \
                class_weight = class_weights, callbacks=[callback, check_point], \
                validation_data=([valid_encode, hla_valid_encode], valid_y_true, w_valid))
        else:
            m1 = pMHCpan.fit(train_encode, train_y_true, \
                epochs=n_epochs, batch_size=batch_size, verbose=2, sample_weight=w_train,\
                class_weight = class_weights, callbacks=[callback, check_point], \
                validation_data=([valid_encode, hla_valid_encode], valid_y_true, w_valid))

        # -----------------------------------------------------------------
        # Prediction performance and Loglikelhood for Train
        # -----------------------------------------------------------------

        if use_nn_encoding:
            y_pred_train = pMHCpan.predict([train_encode, hla_train_encode])
            y_pred_valid = pMHCpan.predict([valid_encode, hla_valid_encode])
        else:
            y_pred_train = pMHCpan.predict(train_encode)
            y_pred_valid = pMHCpan.predict(valid_encode)

        y_pred_train[y_pred_train < MIN_PRED]   = MIN_PRED
        y_pred_train[y_pred_train > 1-MIN_PRED] = 1-MIN_PRED

        y_pred_valid[y_pred_valid < MIN_PRED]   = MIN_PRED
        y_pred_valid[y_pred_valid > 1-MIN_PRED] = 1-MIN_PRED

        df_train['y_pred']  = y_pred_train
        df_valid['y_pred']  = y_pred_valid

        df_train['wt_y_pred']  = y_pred_train
        df_valid['wt_y_pred']  = y_pred_valid

        grp_train = df_train.groupby(['hla', 'peptide', 'y_true'], 
            as_index=False)

        grp_valid = df_valid.groupby(['hla', 'peptide', 'y_true'], 
            as_index=False)

        df_train_max = grp_train['y_pred'].max()
        df_valid_max = grp_valid['y_pred'].max()

        df_train_max.groupby(['y_true'], as_index=False)['y_pred'].describe()
        df_valid_max.groupby(['y_true'], as_index=False)['y_pred'].describe()

        auc_train = roc_auc_score(df_train_max["y_true"], df_train_max["y_pred"])
        auc_valid = roc_auc_score(df_valid_max["y_true"], df_valid_max["y_pred"])

        print('auc with internal encoding one entry per non-binder')
        print([round(auc_train,5), round(auc_valid, 5)])

        ap_train = average_precision_score(df_train_max["y_true"], 
            df_train_max["y_pred"])
        ap_valid = average_precision_score(df_valid_max["y_true"], 
            df_valid_max["y_pred"])

        print('average precision with internal encoding (auc of PR curve)')
        print([round(ap_train,5), round(ap_valid, 5)])

        # -----------------------------------------------------------------
        # multiply predicted density by estiamtes of w
        # -----------------------------------------------------------------

        for ix2 in df_w.index: 
            start_pos1  = df_w.loc[ix2,'start_pos']
            pep_length1 = df_w.loc[ix2,'pep_length']
            w1          = df_w.loc[ix2,'w']
            
            w_pos       = df_train['y_true']==1
            w_pos       = w_pos & (df_train['pep_length']==pep_length1)
            w_pos       = w_pos & (df_train['start_pos']==start_pos1)
            df_train.loc[w_pos,'wt_y_pred'] = df_train.loc[w_pos,'wt_y_pred']*w1
            
            w_pos       = df_valid['y_true']==1
            w_pos       = w_pos & (df_valid['pep_length']==pep_length1)
            w_pos       = w_pos & (df_valid['start_pos']==start_pos1)
            df_valid.loc[w_pos,'wt_y_pred'] = df_valid.loc[w_pos,'wt_y_pred']*w1

        y1_train = (df_train['y_true']==1)
        y1_valid = (df_valid['y_true']==1)

        g1_train = (df_train.loc[y1_train,:]).groupby(['hla','peptide'])
        g1_valid = (df_valid.loc[y1_valid,:]).groupby(['hla','peptide'])

        wt_train = g1_train['wt_y_pred'].transform(lambda x: x / x.sum())
        wt_valid = g1_valid['wt_y_pred'].transform(lambda x: x / x.sum())

        df_train.loc[y1_train, 'weight'] = wt_train
        df_valid.loc[y1_valid, 'weight'] = wt_valid

        # -----------------------------------------------------------------
        # calculate the log-likelihood
        # -----------------------------------------------------------------

        var_grp = ['hla','peptide','y_true']
        df_train_sum = df_train.groupby(var_grp, as_index=False)['wt_y_pred'].sum()
        df_valid_sum = df_valid.groupby(var_grp, as_index=False)['wt_y_pred'].sum()

        y0_train = (df_train_sum['y_true']==0)
        y0_valid = (df_valid_sum['y_true']==0)

        y1_train = (df_train_sum['y_true']==1)
        y1_valid = (df_valid_sum['y_true']==1)

        log_lik_train  = (np.log(binder_weight*df_train_sum.loc[y1_train,'wt_y_pred'])).sum()
        log_lik_train += (np.log(1 - df_train_sum.loc[y0_train,'wt_y_pred'])).sum()

        log_lik_valid  = (np.log(binder_weight*df_valid_sum.loc[y1_valid,'wt_y_pred'])).sum()
        log_lik_valid += (np.log(1 - df_valid_sum.loc[y0_valid,'wt_y_pred'])).sum()

        print('')
        print('log_lik_train and log_lik_valid:')
        print('')
        print(str(log_lik_train)+ '   ' + str(log_lik_valid))
        print('')

        logLiks[idx] = log_lik_train

        # if likelihood decrease, exit the loop. 
        if idx > 0:
            if(log_lik_train + 1e-5*np.abs(log_lik_train) < logLiks[idx-1]):
                warnings.warn('train log likelihood decreased.')
                print('train log decreased')
                print(idx)
                break

        # save the df_w for output at the end of the loop
        # in case the likelihood decreas and we need to remember them
        df_w0  = df_w

        # -----------------------------------------------------------------
        # plot and save losses, and save prediction for validation data
        # -----------------------------------------------------------------

        plot_path  = os.path.join(fig_dir, 
            'track_loss_' + config + '_iter_' + str(idx) + '.png')
        plot_loss(m1, 0, plot_path)

        # -----------------------------------------------------------------
        # Checks for early stopping if the likelihood increases is small
        # for converge_itr consecutative iterations. 
        # ----------------------------------------------------------------- 

        if idx > 0:
            if(log_lik_train < logLiks[idx-1] + converge_e*np.abs(logLiks[idx-1])):
                converge_n = converge_n + 1
            else:
                converge_n = 0
            
            if(converge_n >= converge_itr):
                print('early stopping: convergence of loglikelihood')
                break

    ## end of the EM loop: for idx in range(n_iter):
    end_EM = datetime.now()

    print('log_lik_train across iterations:')
    print(logLiks)
    print('done iteration\n')

    print("Start training =", start_time)
    print("Start EM =", start_EM)
    print("Done EM =", end_EM, "\n")

    fnm_w =  config +  '_w.txt'
    df_w0.to_csv(os.path.join(results_dir, fnm_w), 
            sep='\t', index=False, float_format='%.3e')

    if save_validate_pred: 
        fnm  = config + '_valid_pred_all.txt'
        df_valid.to_csv(os.path.join(results_dir, fnm), 
            sep='\t', index=False, float_format='%.3e')
        subprocess.run(["gzip", os.path.join(results_dir, fnm)])


    print('Finished all iterations of pMHCpan training. Bye!')


# -----------------------------------------------------------------
# parameters
# -----------------------------------------------------------------

parser = argparse.ArgumentParser(
    description='a pan-specific method for HLA-I antigen presentation using deep learning.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "--input-train", 
    type = str,
    dest = "input_train",
    help = "input: path to input file of training data"
)
parser.add_argument(
    "--input-validate", 
    type = str,
    dest = "input_validate",
    help = "input: path to input file of validation data"
)
parser.add_argument(
    "--data-dir", "-D",
    type = str,
    default = "../data/pMHCpan_data",
    help = "directory where data are placed"
)
parser.add_argument(
    "--info-dir", "-I",
    type = str,
    dest = "info_dir",
    default = "../data/NetMHCpan_train",
    help = "directory where HLA allele informaton are placed"
)
parser.add_argument(
    "--fig-dir", "-F",
    type = str,
    dest = "fig_dir", 
    default = "../figures/pMHCpan",
    help = "directory where figures are saved"
)
parser.add_argument(
    "--results-dir", "-R",
    type = str, 
    dest = "results_dir", 
    default = "../results/pMHCpan",
    help = "directory where results are saved"
)
parser.add_argument(
    "--n-layers", "-L",
    type = int,
    default = 1,
    help = "number of layers for encoder and decoder: 1 or 2"
)
parser.add_argument(
    "--hidden-size1", 
    type = int,
    default = 800,
    help = "number of nodes for hidden layer 1"
)
parser.add_argument(
    "--hidden-size2", 
    type = int,
    default = 400,
    help = "number of nodes for hidden layer 2"
)
parser.add_argument(
    "--validating_frac", "-T",
    type = float,
    default = 0.2,
    help = "fraction of data to be used as validation data"
)
parser.add_argument(
    "--batch-size", "-M",
    type = int,
    default = 32,
    help = "batch size used when training"
)
parser.add_argument(
    "--learn-rate",
    type = float,
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
    "--n-epochs", "-e",
    type = int,
    default = 10,
    help = "number of epochs for which to train"
)
parser.add_argument(
    "--n_iter", 
    type = int,
    default = 10,
    help = "number of iterations for which to train"
)
parser.add_argument(
    "--binder_weight",
    type = float,
    default = 5,
    help = "class weight for binders (vs. nonbinder). used only if flag 'use_class_weight'"
)
parser.add_argument(
    "--use_class_weight", 
    dest = "use_class_weight",
    action = "store_true", 
    default = False, 
    help = "Use binder class weights when fitting NN"
)
parser.add_argument(
    "--use_nn_encoding", 
    dest = "use_nn_encoding",
    action = "store_true", 
    default = False, 
    help = "Use a one-layer neural nework to encode HLA and peptide separately"
)
parser.add_argument(
    "--save_validate_pred",
    dest = "save_validate_pred", 
    action = "store_true", 
    default = False, 
    help = "save predictions from validation data"
)



if __name__ == '__main__':
    arguments = parser.parse_args()
    main(**vars(arguments))

