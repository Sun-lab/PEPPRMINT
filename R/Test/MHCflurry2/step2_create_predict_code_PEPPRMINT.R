#---------------------------------------------------------------------------------
# Step 2.Create MHCflurry2.0 code for prediction of PEPPRMINT
# Note: For MA HLA-I EL peptide data 
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Set Parameters ( model configuration)
#--------------------------------------------------------------------------------
hidden.sizes1 = c(200, 400, 800)


#---------------------------------------------------------------------------------
# Create file
#---------------------------------------------------------------------------------
cat("", file="../../../python/run_files/predict_PEPPRMINT_run.sh")

for(split in 0:4){
  cmd0 = sprintf("python3 _prediction_PEPPRMINT.py") 
  for(h in hidden.sizes1){
      cmd = sprintf("%s --input-test-pred ../data/test_data/MHCflurry2/MHCflurry2_test_el_multi_HLA.txt.gz", cmd0)
      cmd = sprintf("%s --test_data_name PEPPRMINT_MHCflurry2", cmd)
      cmd = sprintf("%s --model mixPep_%s_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split%s_Aug5_iter10.h5 ",
                    cmd, h , split)
      cmd = sprintf("%s --m_tag %s_nobal_nocw_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split%s_iter10_Aug5",
                    cmd, h, split)
      cmd = sprintf("%s --results_dir ../results/test_data/MHCflurry2 --input_alist ../data/test_data/MHCflurry2/allelelist_MHCflurry2.txt", 
                    cmd)
      cmd = sprintf("%s --save_all_pred > logfiles/PEPPRMINT_MHCflurry2_%s_bs_32_lr_0.001_e_10_layer_1_dropout_0.5_new_pi_weight_0.5_decriter_2_split%s_Aug5_iter10.log \n",
                    cmd, h, split )
      cat(cmd, file="../../../python/run_files/predict_PEPPRMINT_run.sh", append=TRUE)
  }
  
}

sessionInfo()
q(save="no")

