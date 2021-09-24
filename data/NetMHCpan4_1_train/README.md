# Training Data 

Follow the steps to generate the training data for PEPPRMINT (single allele and multi-allele)

# Step 0: Verify directory set up
If you pull directly from the Github, then the set up should be automatically done. If not, please check that you have the following directories: 
- data/NetMHCpan4_1_train
- data/MA_data
- data/SA_data 
- data/SA_data/by_split

# Step 1: Import data 
Download the training data (NetMHCpan_train.tar.gz) into the 'data/NetMHCpan4_1_train/' directory. Data is provided by the NetMHCpan-4.1 engine located at http://www.cbs.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/. After unzipping the file, there should be five eluted peptide files (c00#_el), where # = 0,...,4. Additionally, there should be an 'allelelist' and 'MHC_pseudo.dat' file, which corresponds to the list of HLA in each sample and the 34 pseudo-sequence for each HLA. 

# Step 2: Generate training data with provided R code 
In the 'R/Train/' directory, there are R scripts to generate training data. 

## step1_extract_SA_data.R
This script will: 
- Extract SA data (eluted peptides that bind to one specific HLA) from each split and put it in the PEPPRMINT training format. Peptides are represented by the 15-length representation (see the README.md on the main page for more details of data format). Duplicated peptides for the same HLA are filtered out. These will be saved in 'data/SA_data/by_split' as 'train_v4_el_single_HLA_#.txt.gz'
- Combine 4 splits to be the training data (train_v4_el_single_HLA.txt.gz) and save in 'data/SA_data'
- Save 1 split to be the validation data (validate_v4_el_single_HLA.txt.gz) and save in 'data/SA_data'
- Reformat 'allelelist' and 'MHC_pseudo.dat' to 'allelelist.txt' and 'MHC_psuedo.txt', respectively. These files are saved to 'data/NetMHCpan4_1_train'.

The '.Rout' file is provided as confirmation of the run script. 
Note: This R script will also create datasets where the SA data is subset to peptides of length 8-11 and 12-15 by split. These will be saved in 'data/SA_data/by_split' with the format 'train_v4_el_8to11_single_HLA_#.txt.gz' or 'train_v4_el_12to15_single_HLA_#.txt.gz'

## step3_extract_MA_data.R
Note: you must have the files provided in 'step1_extract_SA_data.R' in order to run this R script. 

This script will: 
- Extract MA data (eluted peptides that bind to at least one HLA in a sample) from each split and put it in the PEPPRMINT training format. This includes encoding all HLA in the sample for a peptide if it is a binder, and randomly choosing a HLA if it is a non-binder. Peptides are represented by the 15-length representation (see the README.md on the main page for more details of data format). Duplicated peptides for the same sample are filtered out. These will be saved in 'data/MA_data/' with the format 'train_v4_el_multi_only_HLA_#.txt.gz'.
- Add the SA data to MA by split, which is what is used as the training data for PEPPRMINT. These will be saved in 'data/MA_data/' with the format 'train_v4_el_multi_HLA_#.txt.gz'.

The '.Rout' file is provided as confirmation of the run script. 
Note: This R script will also create datasets where the SA+MA data is subset to peptides of length 8-11 and 12-15 by split. These will be saved in 'data/MA_data/' with the format 'train_v4_el_8to11_multi_HLA_#.txt.gz' or 'train_v4_el_12to15_multi_HLA_#.txt.gz'

# Format of Data 
MA data has the following format: 
- 3 columns with headers: peptide, binder, cell_line 
- peptide: the peptide sequence in the 15-length representation. Should be 15 amino acids
- binder : value of 0 or 1 (if a binder = 1, if a non-binder = 0). For the testing set where binding status is unknown, set binder to all 0. 
- cell_line : the sample name. For training data, this is the sample that the peptide is associated with. For testing data, this is the sample to make prediction for. The sample name should match the sample name in the corresponding allele list file, which lists the HLAs in each sample. 

SA data has the following format: 
- 3 columns with headers: peptide, binder, HLA. 
- peptide: the peptide sequence in the 15-length representation. Should be 15 amino acids
- binder : value of 0 or 1 (if a binder = 1, if a non-binder = 0). For the testing set where binding status is unknown, set binder to all 0. 
- HLA : the hla name. For training data, this is the HLA that the peptide is associated with. For testing data, this is the HLA to make prediction for. The HLA name should follow the format "HLA-%##:##", which is 4-digit resolution, % = the letter (A, B, or C), and # = a numerical digit (0 - 9). 

More details are in the README file on the main page.
