#---------------------------------------------------------------------------------
# Step 2.Create MHCflurry2.0 code for prediction of NetMHCpan-4.1
# Note: For MA HLA-I EL peptide data 
#---------------------------------------------------------------------------------
library(stringr)
cell_lines= c("10-002-S1-TISSUE", 
              "11-002-S1-TISSUE",
              "29/14-TISSUE" ,
              "637/13-TISSUE",
              "BCN-018-TISSUE" ,
              "CPH-07-TISSUE" ,
              "CPH-08-TISSUE" ,
              "CPH-09-TISSUE" ,
              "KESKIN_13240-002" ,
              "KESKIN_13240-005" ,
              "KESKIN_13240-006" ,
              "KESKIN_13240-015" ,
              "KESKIN_CP-594_V1" ,
              "KESKIN_DFCI-5283" ,
              "KESKIN_DFCI-5328" ,
              "KESKIN_DFCI-5341" ,
              "KESKIN_H4198_BT187" ,
              "KESKIN_H4512_BT145" ,
              "LEIDEN-004-TISSUE" ,
              "LEIDEN-005-TISSUE" )

allelelist = read.table("../../../data/test_data/MHCflurry2/allelelist_MHCflurry2.txt")
allelelist  = as.data.frame(allelelist )
allelelist$cell_line = str_replace_all(allelelist$V1, "/", "_")
head(allelelist)

#mixPep_100_bs_32_lr_0.001_e_50_layer_1_split0_Oct27.h5

cat("", file="../../../python/run_files/predict_Net4_1_run.sh")

for(c in cell_lines){
  print(c)
  hlas = allelelist$V2[which(allelelist$V1 == c )]
  c_rpl = str_replace_all(c, "/", "_")
  cmd0 = sprintf("../netMHCpan -p MHCflurry2_Net_%s.txt -xls -a", c_rpl) 

  cmd = sprintf("%s %s -xlsfile _predictions/MHCflurry2_Net_%s_pred.xls ", cmd0, 
                toString(hlas), c_rpl)
  cmd = sprintf("%s > _predictions/MHCflurry2_Net_%s_pred.txt \n ", cmd, c_rpl )
      
  cat(cmd, file="../../../python/run_files/predict_Net4_1_run.sh", append=TRUE)
  
}
