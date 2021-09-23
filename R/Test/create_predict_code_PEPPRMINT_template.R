#---------------------------------------------------------------------------------
# Template to create ".sh" file of prediction code for PEPPRMINT
#---------------------------------------------------------------------------------
#Set Parameters 
# string to identify your test set
testFileName = "MHCflurry2" 

# directory of where your test data is located
testFileLoc = "../data/test_data/MHCflurry2/MHCflurry2_test_el_multi_HLA.txt.gz" 

# directory of where your test results should be saved 
testResultsdir = "../results/test_data/MHCflurry2"

# location and name of file that lists all the HLA in each sample 
alistLoc = "../data/test_data/MHCflurry2/allelelist_MHCflurry2.txt"

# set model configuration (PEPPRMINT aggregates acroos 200, 400, 800)
hidden.sizes1 = c(200, 400, 800)

#---------------------------------------------------------------------------------
# Create file
#---------------------------------------------------------------------------------
cat("", file=paste0("../../../python/run_files/", testFileName,
                    "_predict_PEPPRMINT_run.sh"))

for(split in 0:4){
  cmd0 = sprintf("python3 _prediction_PEPPRMINT.py") 
  for(h in hidden.sizes1){
      cmd = sprintf("%s --input-test-pred %s", cmd0, testFileLoc)
      cmd = sprintf("%s --test_data_name PEPPRMINT_MHCflurry2", cmd)
      cmd = sprintf("%s --model mixPep_%s_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split%s_Aug5_iter10.h5 ",
                    cmd, h , split)
      cmd = sprintf("%s --m_tag %s_nobal_nocw_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split%s_iter10_Aug5",
                    cmd, h, split)
      cmd = sprintf("%s --results_dir %s --input_alist %s", 
                    cmd, testResultsdir, alistLoc)
      cmd = sprintf("%s --save_all_pred > logfiles/PEPPRMINT_%s_%s_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split%s_Aug5_iter10.log \n",
                    cmd, testFileName, split )
      cat(cmd, file=paste0("../../../python/run_files/", testFileName,
                           "_predict_PEPPRMINT_run.sh"), append=TRUE)
  }
  
}

sessionInfo()
q(save="no")

