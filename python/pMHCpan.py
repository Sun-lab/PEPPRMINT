
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

# input_train  = "train_el_single_HLA_9aa_0.txt"
# input_test   = "test_el_single_HLA_9aa_0.txt"
# data_dir     = "../data/"
# fig_dir      = "../figures"
# results_dir  = "../results"

# testing_frac = 0.2

# hidden_size1 = 400
# hidden_size2 = 200
# n_layers     = 1

# n_epochs     = 100
# batch_size   = 32
# learn_rate   = 1e-3

def main(input_train, input_test, olabel,
    data_dir, info_dir, fig_dir, results_dir, 
    testing_frac, n_layers,  hidden_size1, hidden_size2,  
    n_epochs, batch_size, learn_rate,  binder_weight=21,
    use_class_weight = True, save_model=True, save_test_pred = True):
    
    start = datetime.now()
    print("start time :", start)

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
        config = config + '_used_class_weight'

    config = config + '_' + olabel

    config

    # -----------------------------------------------------------------
    # read in data
    # -----------------------------------------------------------------

    file_name = os.path.join(data_dir, input_train)
        
    print('input file:')
    print(file_name)
    print('')
        
    df = pd.read_csv(file_name, sep='\t', header=0)
    df.shape
    df.iloc[0:3,0:5]

    # -----------------------------------------------------------------
    # split_training_testing_data
    # -----------------------------------------------------------------

    if input_test:
        trainX    = df
        file_test = os.path.join(data_dir, input_test)
        testX     = pd.read_csv(file_test, sep='\t', header=0)
        print('file name for testing data is provided, will not split data\n')
    else:
        print('file name for testing data is not provided, will split data\n')
        dfT = df.T
        trainX, testX = train_test_split(dfT, test_size=testing_frac, random_state=1999)

    print('training and testing data dimension:')
    print(trainX.shape)
    print(testX.shape)

    print('trainX[0:2,]:')
    print(trainX.iloc[0:2,])
    print('testX[0:2,]:')
    print(testX.iloc[0:2,])

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
    #print(hla_encoding['BoLA-100901'])

    # -----------------------------------------------------------------
    # encoding peptide and hla data
    # -----------------------------------------------------------------

    train_encode = list()
    test_encode  = list()

    for i in trainX.index:
        pep_encode = encode_peptide(trainX.iloc[i,0])
        hla_encode = hla_encoding[trainX.iloc[i,2]]
        train_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())

    for i in testX.index:
        pep_encode = encode_peptide(testX.iloc[i,0])
        hla_encode = hla_encoding[testX.iloc[i,2]]
        test_encode.append(np.concatenate((pep_encode,hla_encode)).flatten())

    train_encode = np.array(train_encode)
    test_encode  = np.array(test_encode)

    print('train_encode[0:2,]:')
    print(train_encode[0:2,])
    print(train_encode.shape)

    print('test_encode[0:2,]:')
    print(test_encode[0:2,])
    print(test_encode.shape)
    print('')

    # -----------------------------------------------------------------
    # model setup
    # -----------------------------------------------------------------
    #define class weights    
    bind_class_weight = {
        0: 1.,
        1: binder_weight}

    if use_class_weight:
        print('')
        print("----------------------------------------------------")
        print("~~~~~ binder class weight used in NN ")
        print("----------------------------------------------------")
        print(bind_class_weight)
        print('')

    input1 = keras.Input(shape=(train_encode.shape[1],), name='input')
    hidden = layers.Dense(hidden_sizes[0], activation='relu', name='hidden0')(input1)
    hidden = layers.Dropout(rate=0.5, name='dropout_hidden0')(hidden)

    if n_layers > 1: 
        hidden = layers.Dense(hidden_sizes[1], activation='relu', name='hidden1')(hidden)
        hidden = layers.Dropout(rate=0.5, name='dropout_hidden1')(hidden)

    prob = layers.Dense(1, activation='sigmoid', name='output')(hidden)

    pMHCpan = keras.Model(input1, prob)
    pMHCpan.summary()

    # -----------------------------------------------------------------
    # model fitting
    # -----------------------------------------------------------------

    adam1 = keras.optimizers.Adam(lr=learn_rate)    
    pMHCpan.compile(optimizer=adam1, loss='binary_crossentropy', metrics=['accuracy'])

    train_y_true = np.array(trainX['binder'])
    test_y_true  = np.array(testX['binder'])

    if use_class_weight:
        print('')
        print("~~~~~ fit NN with class weights")
        m1 = pMHCpan.fit(train_encode, train_y_true, epochs=n_epochs, 
        batch_size=batch_size, verbose=2, class_weight = bind_class_weight,
        validation_data=(test_encode, test_y_true))
    else: 
        print('')
        print("~~~~ fit NN without class weights")
        m1 = pMHCpan.fit(train_encode, train_y_true, epochs=n_epochs, 
        batch_size=batch_size, verbose=2, 
        validation_data=(test_encode, test_y_true))

    # -----------------------------------------------------------------
    # Prediction performance and Loglikelhood for Train
    # -----------------------------------------------------------------
    train_pred = pMHCpan.predict(train_encode)
    out_train = pd.DataFrame(data={'y_pred': train_pred.flatten(), 
            'y_true': trainX['binder']})

    #loglikelhiood  -- if would like to calculate

    #y0_train = (out_train['y_true']==0)
    #y1_train = (out_train['y_true']==1)

    #if use_class_weight:
    #    log_lik_train = (np.log(binder_weight* out_train.loc[y1_train,'y_pred'])).sum()
    #else:
    #    log_lik_train = (np.log(out_train.loc[y1_train,'y_pred'])).sum()
            
    #log_lik_train = log_lik_train + (np.log(1 - out_train.loc[y0_train,'y_pred'])).sum()
        
    #print('')
    #print('~~~~~~ log_lik_train ~~~~~~')
    #print('')
    #print(str(log_lik_train))

    # Performance of TRaining
    auc_pMHC_train   = roc_auc_score(out_train["y_true"], out_train["y_pred"])
    pr_pMHC_train    = average_precision_score(out_train["y_true"], out_train["y_pred"])

    print('')
    print("----------------------------------------------------")
    print('TRAIN: Evaluation')
    print('auc ROC, auc PR')
    print("----------------------------------------------------")
    print([round(auc_pMHC_train,3), round(pr_pMHC_train,3)])
    print('')

    # -----------------------------------------------------------------
    # plot and save losses, and save prediction for testing data
    # -----------------------------------------------------------------

    plot_path  = os.path.join(fig_dir, 'track_loss_' + config + '.png')
    plot_loss(m1, 1, plot_path)

    loss = pd.DataFrame.from_dict(m1.history)
    fnm  =  config + '_loss.txt'
    loss.to_csv(os.path.join(results_dir, fnm), sep='\t', 
        index=False, float_format='%.3e')

    if save_test_pred:
        test_pred = pMHCpan.predict(test_encode)
        out_test = pd.DataFrame(data={'pred': test_pred.flatten(), 
            'true': testX['binder']})

        fnm      = config + '_test_pred.txt'
        out_test.to_csv(os.path.join(results_dir, fnm), sep='\t', 
            index=False, float_format='%.3e')
        subprocess.run(["gzip", os.path.join(results_dir, fnm)])

        auc_pMHC   = roc_auc_score(out_test["true"], out_test["pred"])
        pr_pMHC    = average_precision_score(out_test["true"], out_test["pred"])

        print("----------------------------------------------------")
        print('VALID: Evaluation')
        print('auc ROC, auc PR')
        print("----------------------------------------------------")
        print([round(auc_pMHC,3), round(pr_pMHC,3)])

    

    # -----------------------------------------------------------------
    # save the model
    # -----------------------------------------------------------------

    if save_model:
        model_path  = os.path.join(results_dir, config + '.h5')
        pMHCpan.save(model_path)

    end = datetime.now()

    print("----------------------------------------------------")
    print('Done training pMHCpan!')
    print("----------------------------------------------------")
    print("Started training =", start)
    print("Finished training =", end, "\n")
# -----------------------------------------------------------------
# parameters
# -----------------------------------------------------------------

parser = argparse.ArgumentParser(
    description='a pan-specific method for antigen presentation using deep learning.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "--input-train", 
    type = str,
    dest = "input_train",
    help = "input: path to input file of training data"
)
parser.add_argument(
    "--input-test", 
    type = str,
    dest = "input_test",
    help = "input: path to input file of testing data"
)
parser.add_argument(
    "--data-dir", "-D",
    type = str,
    default = "../data/SA_data",
    help = "directory where data are placed"
)
parser.add_argument(
    "--info-dir", "-I",
    type = str,
    dest = "info_dir",
    default = "../data/NetMHCpan4_1_train",
    help = "directory where HLA allele informaton are placed"
)
parser.add_argument(
    "--fig-dir", "-F",
    type = str,
    default = "../figures/SA",
    help = "directory where figures are saved"
)
parser.add_argument(
    "--results-dir", "-R",
    type = str,
    default = "../results/SA",
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
    "--testing-frac", "-T",
    type = float,
    default = 0.2,
    help = "fraction of data to be used as testing data"
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
    default = 50,
    help = "number of epochs for which to train"
)
parser.add_argument(
    "--binder_weight",
    type = float,
    default = 21,
    help = "class weight for binders (vs. nonbinder). used only if flag 'use_class_weight'"
)
parser.add_argument(
    "--use_class_weight", 
    dest = "use_class_weight",
    action = "store_true", 
    help = "Use binder class weights when fitting NN"
)
parser.add_argument(
    "--save_model",
    dest = "save_model",
    action = "store_true",
    help = "save the model to a h5 file"
)
parser.add_argument(
    "--save_test_pred",
    dest = "save_test_pred",
    action = "store_true",
    help = "save validation predictions"
)
parser.set_defaults(save_model = True, save_test_pred = True)


if __name__ == '__main__':
    arguments = parser.parse_args()
    main(**vars(arguments))

