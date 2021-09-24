# MA Test Data (provided by MHCflurry-2.0)

We share the MA test data (Data_S1.csv.gz) complied by MHCflurry-2.0 under the CC BY 4.0 license and provided at https://data.mendeley.com/datasets/zx3kjzc3yx/3. 

# R/Test/MHCflurry2/step0_create_MA_test_data.R
Generates the MA test data in the format to be used for PEPPRMINT. This script does the following: 
- Formats MA test data. Each peptide is encoded for each HLA in the sample and represnted by the 15-length sequence. These are saved in 'data/test_data/MHCflurry2/' as 'MHCflurry2_test_el_multi_HLA.txt'
- Generates an allele list for the test data, saved in 'data/test_data/MHCflurry2/' as 'allelelist_MHCflurry2.txt'

The '.Rout' file is may be used for reference or confirmation. 

Note: for a general template of formating test data for PEPPRMINT, please see 'R/Test/PEPPRMINT_pred_performance_template.R'

# Predict peptide presentation for MA test set 
To make prediction using PEPPRMINT, please see "predict_PEPPRMINT_run.sh" located in 'python/run_files/' directory. 

To make prediction using NetMHCpan-4.1, please see "predict_Net4_1_run.sh" located in 'python/run_files/' directory. 

# Format 
MA data has the following format: 
- 3 columns with headers: peptide, binder, cell_line 
- peptide: the peptide sequence in the 15-length representation. Should be 15 amino acids
- binder : value of 0 or 1 (if a binder = 1, if a non-binder = 0). For the testing set where the true binding status is unknown, set binder to all 0. 
- cell_line : the sample to make prediction for. The sample name should match the sample name in the corresponding allele list file, which lists the HLAs in each sample. 


  
