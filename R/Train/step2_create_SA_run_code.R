#---------------------------------------------------------------------------------
# Step 2. Create Terminal code to run pMHCpan
# Note: For SA HLA-I EL peptide data 
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Set parameters
#---------------------------------------------------------------------------------
hidden.sizes1 = c(100, 200, 400, 800)
Date = "Jul26"

cat("", file="../../python/run_files/SA_training_v4_run.sh")

# SA train data
i1 = "train_v4_el_single_HLA.txt.gz"
# SA validation data
i2 = "validate_v4_el_single_HLA.txt.gz"

#---------------------------------------------------------------------------------
# Example to create run code 
# Model training here is not balanced data, including duplicates, 
#    and not using class weight option
#---------------------------------------------------------------------------------
cmd0 = sprintf("python3 pMHCpan.py --input-train %s --input-test %s", i1, i2)
  
for(l in 1:1){
  for(h1 in hidden.sizes1){
    h2 = h1/2
      
    if(l==1){
      config = sprintf("%d_layer_%d", h1, l)
    }else if(l==2){
      config = sprintf("%d_%d_layer_%d", h1, h2, l)
    }
      
  cmd = sprintf("%s --hidden-size1 %d --hidden-size2 %d", cmd0, h1, h2)
  cmd = sprintf("%s -L %d --olabel %s  --save_test_pred --save_model > logfiles/pMHCpan_v4_%s_nobal_nocw_full_%s.log \n", 
              cmd, l, paste("full", Date, sep = "_"), config, Date)
  cat(cmd, file="../../python/run_files/SA_training_v4_run.sh", append=TRUE)
  }
}




sessionInfo()
q(save="no")

#---------------------------------------------------------------------------------
# Example code using cw
#---------------------------------------------------------------------------------

# python3 pMHCpan.py --input-train train_v4_el_single_HLA.txt.gz --input-test train_v4_el_single_HLA_4.txt.gz 
#   --hidden-size1 800 --hidden-size2 400 -L 1 --olabel full_Jul26 
#   --binder_weight 18 --use_class_weight --save_test_pred --save_model > logfiles/pMHCpan_v4_800_layer_1_classweight18_full_Jul26.log
