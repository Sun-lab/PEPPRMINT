# Results 

See below for description of directories: 

# PEPPRMINT directory
Trained MA models (.h5) and corresponding files (_pi.txt) are saved here. These models are used to make prediction on testing data 

#SA directory
Trained SA model used to initialize the MA models 

# test_data directory
Test set predictions using PEPPRMINT or NetMHCpan-4.1 are saved here. Files that end in 'test_pred_all.txt.gz' are all the predictions, i.e. each peptide predicted presentation by each HLA in the sample. Files that end in 'test_pred.txt' are the maximum predicted value for each peptide (across all HLA in the sample). 

## Prediction files using '_prediction_PEPPRMINT.py'
A prediction file using '_prediction_PEPPRMINT.py' has 4 colums (y_pred_mix,	y_true,	peptide,	hla). 
- "y_pred_mix" is the prediction for the peptide and HLA using the MA trained model (which is used for PEPPRMINT)
- "y_true" is the true binding status, or if unknown, what was given as input (we recommnded to set all to "0")
- "hla" is the HLA in the sample that the prediction is for 
- "peptide" is the concatenated name of the sample and the peptide that the prediction is for in the 15-length representation
- "hla" is the HLA in the sample that the prediction is for

A snapshot of the prediction file is below: 
y_pred_mix	y_true	peptide	hla
6.111e-01	1	10-002-S1-TISSUE;QLEDXXXEXXXALKY	HLA-A02:01
6.480e-01	1	10-002-S1-TISSUE;QLEDXXXEXXXALKY	HLA-A31:01
9.016e-01	1	10-002-S1-TISSUE;QLEDXXXEXXXALKY	HLA-B13:02
5.020e-02	1	10-002-S1-TISSUE;QLEDXXXEXXXALKY	HLA-B58:01
7.419e-03	1	10-002-S1-TISSUE;QLEDXXXEXXXALKY	HLA-C06:02

## Prediction files using a test dataset in the Training 'PEPPERMINT.py'
If a test set is specified using "--input_test" in 'PEPPRMINT.py', then the resulting prediction file has 6 columns(sample,	hla,	peptide,	y_true,	y_pred_pMHC,	y_pred). Note dataset specified for "--input_test" is not used in any way during the training. 
- "sample" is the sample name. This is the HLA for SA test data
- "hla" is the HLA in the sample that the prediction is for 
- "peptide" is the peptide that the prediction is for in the 15-length representation
- "y_true" is the true binding status, or if unknown, what was given as input (we recommended to set all to "0")
- "y_pred_pMHC" is the prediction for the peptide and HLA using the SA trained model 
- "y_pred" is the prediction for the peptide and HLA using the MA trained model (which is used for PEPPRMINT)

A snapshot of the data is as follows: 
sample	hla	peptide	y_true	y_pred_pMHC	y_pred
HLA-A02:02	HLA-A02:02	AAADXXXIXXXVNFL	1	7.326e-01	8.956e-01
HLA-A02:02	HLA-A02:02	AAFEXXXFDXXIHQV	1	9.063e-01	6.995e-01
HLA-A02:02	HLA-A02:02	ADASXXXLXXXLKKV	1	1.497e-01	5.794e-01
HLA-A02:02	HLA-A02:02	AEISXXQIHQXSVTD	1	3.370e-02	3.441e-02




  
