# PEPPRMINT: PEPtide PResentation using a MIxture model and Neural neTwork

PEPPRMINT was designed to exploit the newly available mass spectrometry data as training data to predict peptide presentation by HLAs. In multi-allele mass spectrum data, the peptides can be presented by any one of multiple HLA alleles. PEPPRMINT models these peptides by mixture model where the mixture componenets are neural networks with the same architecure but different inputs. Other than NetMHCpan-4.1, it is the only neural network method that uses multi-allele mass spectrum data for training. 

Next we describe two situations to use our software. One is to make prediction using our pre-trained neural networks. The other is to train the neural networks using another training data.

# Prediction 

The final prediction score for PEPPRMINT is an average of 15 models (3 neural network configurations and 5 splits of training data). Currently, users need to (1) run the prediction of the 15 models in python and then (2) aggregate them using the provided R code to obtain the final prediction. 

## Step 1. Run prediction

The code to make prediction is ```_prediction_PEPPRMINT.py```, which reads input from a tab-delimitated text file with at least three columns. The first row of this file is treated as headers, though the code identifies each column by its order rather than the header. The first three columns should be prepared as follows. 

- The first column is the peptide sequence of 15 amino acids long. 

- The second column is either binding indicator or label of the corresponding somatic mutation. 
	- When evaluating the performance of PEEPRMINT using a dataset with binder status, it is 1 if binder, and 0 if non-binder. If the binder information is unknown, the second column can be set to be 0 for all entries. 
	- When using the code to prioritize neoantigens for cancer vaccine, each somatic mutation corresponds to multiple peptides, and we take the maximum across the peptides for each somatic mutation. In such case, the second column is a label for the somatic mutation. In this  case, the option ```--neoantigen``` should be used when running the code. 
	
- The third column is the corresponding sample of the peptide. The name of the sample will be used to identify the corresponding HLA alleles from another file specified by the ```input_alist``` parameter. 

See the files in folder 'data/_toy_data' folder some examples of input files. Another example of input file for neoantigen analysis can be found in folder 'data/Riaz_2017/riaz_peptide_mut.txt'

Below is the sample python code to run prediction for one of the 15 models. Note "TESTNAME" is a one word indicator to name your testing set that does not include any "_".  

```{python}
python3 _prediction_PEPPRMINT.py --input_test_pred ../data/test_data/your_test_file.txt --test_data_name TESTNAME --model MA_200_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split0_Aug5_iter10.h5  --m_tag 200_split0 --results_dir ../results/PEPPRMINT --input_alist ../data/test_data/MHCflurry2/your_allelelist_file.txt --multi_allele --save_all_pred > logfiles/logfile_name.log 
```

Repeat the above code, except change the "--model"  for the following 15 model configurations where 'xxx' stands for `_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_'. Also set the m_tag according to the input model. For example, m_tag 200_split0 corresponds to the model with 200 nodes in the hidden layer and trained on the split 0 of the training data. 

* MA_200_xxx_split0_Aug5_iter10.h5
* MA_400_xxx_split0_Aug5_iter10.h5 
* MA_800_xxx_split0_Aug5_iter10.h5  
* MA_200_xxx_split1_Aug5_iter10.h5 
* MA_400_xxx_split1_Aug5_iter10.h5  
* MA_800_xxx_split1_Aug5_iter10.h5  
* MA_200_xxx_split2_Aug5_iter10.h5  
* MA_400_xxx_split2_Aug5_iter10.h5 
* MA_800_xxx_split2_Aug5_iter10.h5  
* MA_200_xxx_split3_Aug5_iter10.h5 
* MA_400_xxx_split3_Aug5_iter10.h5  
* MA_800_xxx_split3_Aug5_iter10.h5 
* MA_200_xxx_split4_Aug5_iter10.h5 
* MA_400_xxx_split4_Aug5_iter10.h5 
* MA_800_xxx_split4_Aug5_iter10.h5

The output will be a '*_test_pred_all.txt.gz' file with a prediction for each each peptide and cell line (i.e. across all possible HLA in the cell line for HLA-I). The columns of the output file is as follows: 

- peptide : Name of cell line and the peptide 
- hla : HLA in cell line 
- y_true : true binding of the peptide. (Note, if do not have the true binding status, y_true = 0 for all peptides)
- y_pred_mix : prediction of binding probability for the peptide and HLA

## Step 2. Aggregate Predictions for PEPPRMINT

In the 'PEPPRMINT_pred_performance_template.R' code (located in R/Test/), the maximum prediction for each peptide across all possible HLA is found. These maximums are then aggregated across the 15 models through an average and saved as 'cv_score'. These aggregated scores are the final PEPPRMINT prediction. 

## Parameter specification in PEPPRMINT.py

- input_test_pred : test data file name. File should be a '.txt' file and delimited by '\t'. 
- input_alist : name of file with a list of the HLA-I alleles for all the samples or cell lines. Default is 'allelelist.txt' for HLA-I
- input_MHCseq : name of file with the 34 pseudo-sequence for all HLA alleles. Default is 'MHC_psuedo.txt' for HLA-I
- test_data_name : label for the test data to use in the output
- model : specification of '.h5' for the model to use for prediction
- m_tag: label for the model used for prediction 
- data_dir (-D): pathway to directory where datasets are located. Default is '../data/test_data'.
- fig_dir (-F): pathway to directory where output figures will be saved. Default is '../figures/test_data'
- results_dir (-R): pathway to directory where results will be saved (including model if 'save_model = True'). Default is '../results/test_data'
- info_dir (-I): pathway to directory where the allele list and pseudo-sequence files are located. Default is '../data/NetMHCpan4_1_train'
- save_all_pred : option to save all predictions for all combinations of peptide, HLA, and binding core (if applicable) in cell line. File will be saved with suffix 'test_pred_all'. Default is 'False'.

## Supporting file format
- The provided default file for the 34 position pseudo sequence is 'MHC_psuedo.txt' and additional HLA can be added. If want to use a different file, it must be delimited by ' ' and have no header. The two columns are the HLA name and the 34 pseudo sequence.  

- The provided default file for that lists the possible HLA-I per sample is 'allelelist.txt'. New files should be '.txt' file delimited by ' ' (space) and have no header. The two columns are the cell line name and possible HLA in the cell line. The list of HLA are separated by commas. The HLA-I should be 4 resolution and follow the the following nomenclature "HLA-\$##:##", where \$ = (A, B, or C) and # = integer from 0 to 9.


# Training Multi-allele models

The neural networks in PEPPRMINT were trained with the same data as NetMHCpan-4.1. The reference for the paper is as follows: Alvarez, B., Reynisson, B., Barra, C., Buus, S., Ternette, N., Connelley, T., ... & Nielsen, M. (2019). NNAlign_MA; MHC peptidome deconvolution for accurate mhc binding motif characterization and improved t-cell epitope predictions. Molecular & Cellular Proteomics, 18(12), 2459-2477.

## Train PEPPRMINT.py
Sample code to train a one layer model with 100 hidden nodes and a given test validation set is below. Model will be saved as a '.h5' file. The training is initlized by previously trained SA HLA-I, respectively. An optional file with maximum prediction for a peptide the given validation set across HLA in the cell line can be saved in a '*_test_pred.txt' using the 'val_pred_data' option. Note: the version of Tensorflow used the train this model must match the version of Tensorflow that was used to train the initializing model. The following is the python call code to train a HLA-I MA model. 

```{python}
python3 PEPPRMINT.py --input_train train_MA_sample_data.txt --input_valid valid_MA_sample_data.txt --input_test test_MA_sample_data.txt --init_model ../results/SA/pMHCpan_800_bs_32_lr_0.001_e_50_layer_1_full_Jul26.h5 --hidden-size1 800 --hidden-size2 400 -L 1 --olabel split4_Aug5  --n_epochs 10 --decr_iter 2 --save_model --save_all_pred > logfiles/sample_800_split4_Aug5.log 
```

## Parameter specification for PEPPRMINT.py
The parameters for the neural networks of PEPPRMINT are as follows. 
- input_train: specify the training dataset following format described below
- input_valid : specify the validation test dataset following same format as training dataset
- input_test : specify an indpendent test dataset that will use external encoding. Note, this dataset is not used in the training!
- input_alist : list of HLA alleles for all the samples. Default file is 'allelelist.txt' for HLA-I. The format of the file is specified in 'Data Format' section. 
- input_MHCseq : 34 position pseudo-sequence of different HLA-I alleles. Default file is 'MHC_pseudo.txt' for HLA-I. See above for formating.
- init_model : specification of '.h5' file with SA HLA-I model weights which will initialize the training. Default is 'pMHCpan_800_bs_32_lr_0.001_e_50_layer_1_full_Jul26.h5' (though don't have to include a defualt)
- data_dir (-D): path to directory where datasets are located. Default is '../data/MA_data'.
- info_dir (-I): path to directory where allele list and HLA pseudo-sequence files are located. Default is "../data/NetMHCpan4_1_train"
- fig_dir (-F): path to directory where output figures will be saved. Default is '../figures/PEPPRMINT'
- results_dir (-R): path to directory where results will be saved (including model if 'save_model = True'). Default is '../results/PEPPRMINT' 
- n_layers (-L): specify number of layers for encodera nd decoder as 1 or 2 in model structure. Default is 1. 
- hidden_size1 : specify number of nodes in the first layer. Default is 800. 
- hidden_size2: if applicable, specify the number of nodes in the second layer. Default is 400. Even if specified, this will not be used if n_layers = 1.  
- n_epochs (-e): number of training epochs. Default is 50 
- batch_size (-M) : batch size default is 32 
- learn_rate: learning rate default is 0.001 
- olabel: optional string to add to all saved output names. 
- n_iter : Maximim number of iterations of the EM algorithm. Default is 15. 
- drop_rate : Drop out rat used in the NN layers. Default is 0.5
- decr_iter : number of times the log likelihood can decrease before stopping the training early
- use_class_weight: option to use binder class weights when fitting NN. Default is False
- converge_e : specified amount considered as converging of likelihood. Default is 
- converge_itr : maximum number of iterations where log likelihood has not increased more than 'converge_e' before early stopping of training. 
- save_model : option to save the final model to a '.h5' file. Use '--save_model' option if want to save model. 
- save_all_iterations: save all the iteration results
- save_all_pred: option to save all predictions in validation set (instead of just maximum). Use '--save_all_pred' option in command line. 

# Training Single allele models (SA)

The single allele (SA) method trained on HLA-I mass spectrum that is used to initialize our multi-allele method PEPPRMINT. The default model structure is one dense layer, though there is the option to train a two-layer model. The model is trained with a binary cross-entropy loss function and ADAM optimizer. Below are the details to train these SA models. 

## Parameter specification
The parameters include the following and are defined the same way as the MA methods PEPPRMINT unless otherwise noted.
- input_train, data_dir (-D), fig_dir (-F), results_dir (-R),  n_layers (-L), hidden_size1, hidden_size2, n_epochs (-e), batch_size (-M), learn_rate, olabel, save_model, val_data_pred 
- input_test : specify the validation test dataset following same format as training dataset. Could also specify a 'testing_frac' as the validation test set. 
- testing_frac (-T): fraction of training dataset that will be used as a validation test set. Is not needed if specified 'input_test'. Default is 0.2. 

## Data format
- HLA-I input data files for training and testing: '.txt' files delimited by '\t' and have no header. The columns of the datasets are peptide sequence in the 15-length representation, binder indicator (=1 if binder, 0 if non-binder), and specific HLA that the peptide binds to. The HLA should be 4 resolution and follow the the following nomenclature "HLA-\$##:##", where \$ = (A, B, or C) and # = integer from 0 to 9. See 'pMHCpan_sample_data.txt' in the 'data' folder for a example of the data format. 

- Format of HLA pseudo-sequence or cell line list files: please refer to the 'Data Format' section under the multi-allele training section above. 

## Train Model 
Sample code to train a one layer model with 800 hidden nodes and a given validation set is below. Model will be saved as '.h5' file. Predictions for the given validation set will be saved in a '*_test_pred.txt' file. 

```{python}
python3 pMHCpan.py --input-train train_ds.txt --input-test test_ds.txt -L 1 --hidden-size1 100 --hidden-size2 50 --output_tag sampleRun --val_data_pred --save-model > logfiles/pMHCpan_sampleRunlog
```
  
