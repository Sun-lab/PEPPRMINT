# SA Test Data (provided by NetMHCpan-4.1)

Follow the steps below to generate the SA test data provided by NetMHCpan-4.1

# Step 1: Download SA test data 
Download the testing data into the 'data/test_data/NetMHCpan4_1' directory. Data is provided by the NetMHCpan-4.1 engine located in the "Evaluation data, MS Ligands" section on http://www.cbs.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/. These test data are provided as 20 separate files, one for each HLA (i.e. HLA-A02:02, ..., HLA-C12:03). 

# Step 2: Check name and type of file
Downloaded files should be '.txt' files and the ":" replaced with "_"(HLA-A02:02 is saved as HLA-A02_02.txt). These can be done in the Finder window.

# Step 3: Run R/Test/NetMHCpan4_1/step0_combine_SA_test_set.R
Generates the SA test data in the format to be used for PEPPRMINT. This script does the following: 
- Formats SA test data. Each peptide is encoded for each HLA in the sample and represnted by the 15-length sequence. These are saved in 'data/test_data/NetMHCpan4_1/' as ' 'NetMHCpan41_test_el_single_HLA.txt'

The '.Rout' file is may be used for reference or confirmation. 

Note: for a general template of formating test data for PEPPRMINT, please see 'R/Test/PEPPRMINT_pred_performance_template.R'

# Step 4: Peptide Prediction for SA test set 
In PEPPRMINT.py, there is an option to include an independent test set, where the peptide presentation is predicted at the end of training. This test set is not used in the training. We predicted SA test set using this method. To make predictions, code similar to below should be used. See the 'README.md' on the main page for more details. 

python3 _prediction_PEPPRMINT.py --input-test-pred ../data/test_data/NetMHCpan4_1/NetMHCpan41_test_el_single_HLA.txt.gz --test_data_name PEPPRMINT_NetMHCpan41 --model MA_200_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split0_Aug5_iter10.h5  --m_tag 200_split0 --results_dir ../results/PEPPRMINT --save_all_pred > logfiles/logfile_name.log 

For NetMHCpan-4.1 prediction, we used the portable version of NetMHCpan-4.1 located at http://www.cbs.dtu.dk/services/NetMHCpan/. 

# Format of Data
SA data for PEPPRMINT has the following format: 
- 3 columns with headers: peptide, binder, HLA. 
- peptide: the peptide sequence in the 15-length representation. Should be 15 amino acids
- binder : value of 0 or 1 (if a binder = 1, if a non-binder = 0). For the testing set where binding status is unknown, set binder to all 0. 
- HLA : the HLA name to make prediction for. The HLA name should follow the format "HLA-%##:##", which is 4-digit resolution, % = the letter (A, B, or C), and # = a numerical digit (0 - 9). 


  
