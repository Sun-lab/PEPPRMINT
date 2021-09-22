# SA Data 
 
This is the folder that will have all the single-allele mass spectrum (MA) training data generated for the burn-in period used to intitialize the neural networks in PEPPRMINT. 

To generate SA data, please follow the steps outlined in 'README_train_data.md' in the 'data/NetMHCpan4_1_train' directory. 

SA data has the following format: 
- 3 columns with headers: peptide, binder, HLA. 
- peptide: the peptide sequence in the 15-length representation. Should be 15 amino acids
- binder : value of 0 or 1 (if a binder = 1, if a non-binder = 0). For the testing set where binding status is unknown, set binder to all 0. 
- HLA : the hla name. For training data, this is the HLA that the peptide is associated with. For testing data, this is the HLA to make prediction for. The HLA name should follow the format "HLA-%##:##", which is 4-digit resolution, % = the letter (A, B, or C), and # = a numerical digit (0 - 9). 
  
