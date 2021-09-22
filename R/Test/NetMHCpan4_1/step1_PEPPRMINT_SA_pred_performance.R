# ---------------------------------------------------------------------------
# Analyze PEPPRMINT prediction results 
# NetMHCpan-4.1 test data 
# ---------------------------------------------------------------------------
library(data.table)
library(PRROC)
library(stringr)

# ---------------------------------------------------------------------------
# Import data 
# ---------------------------------------------------------------------------
dir0 = "../../../results/test_data/NetMHCpan4_1/"
fall = list.files(dir0, 
                  pattern = "test_pred_all.txt.gz")
fall

# ---------------------------------------------------------------------------
#extract specific HLA
# ---------------------------------------------------------------------------
p =  fread(paste0(dir0, fall[1]))

hits_tab  = table( p$sample, p$y_true)
hits_tab

hits_tab1 = hits_tab[,2]
hits_tab1


# ---------------------------------------------------------------------------
#function for each iteration
# ---------------------------------------------------------------------------
SAel_all_max = NULL 

for(f in fall){
  print("*-------------------------------------------------------------------------*")
  it = str_split(f, "_")[[1]][19]
  #print(it)
  
  p_pred = fread(paste0(dir0, f))
  
  print(dim(p_pred))
  head(p_pred)
  
  #Remove duplicated peptide within the same HLA
  p_pred$pep_line = paste(p_pred$sample, p_pred$peptide, sep= ";")
  p_pred = p_pred[!duplicated(p_pred$pep_line)]
  
  # ---------------------------------------------------------------------------
  # Overall split/config performance (validate python results)
  # ---------------------------------------------------------------------------
  # calculate ROC curve 
  roc1 = roc.curve(scores.class0 = p_pred$y_pred, weights.class0 = p_pred$y_true, 
                   curve=TRUE)
  roc1
  #plot(roc1)
  
  #pr1  = pr.curve(scores.class0 = d50$y_pred_mix, weights.class0 = d50$y_true,
  #                curve=TRUE)
  #pr1
  #plot(pr1)
  
  # PPV
  temp1 <- p_pred[with(p_pred,order(-y_pred)),]
  hits = table(p_pred$y_true)[2][[1]]
  temp2 = temp1[1:hits,]
  dim(temp2)
  
  ppv = sum(temp2$y_true)/hits
  ppv
  
  print(paste0("AUC, PPV for config: ", f))
  print(paste("AUC =", round(roc1$auc,3)," PPV = ", round(ppv,3), sep = " "))
  
  # ---------------------------------------------------------------------------
  # Combine all split/configs to one dataset
  # ---------------------------------------------------------------------------
  colnames(p_pred)[which(colnames(p_pred)=="y_pred")] = paste("y_pred_MA",
                                                            str_split(f, "_")[[1]][2],
                                                            str_split(f, "_")[[1]][19], 
                                                            sep = "_")
  colnames(p_pred)[which(colnames(p_pred)=="y_pred_pMHC")] = paste("y_pred_pMHC",
                                                              str_split(f, "_")[[1]][2],
                                                              str_split(f, "_")[[1]][19], 
                                                              sep = "_")
  
  if(str_split(f, "_")[[1]][2]=="200" & str_split(f, "_")[[1]][19]=="split0"){
    SAel_all_max = p_pred
    print(dim(SAel_all_max))
  }else{
    SAel_all_max = merge(SAel_all_max, p_pred, 
                    by = intersect(colnames(SAel_all_max), colnames(p_pred)))
    print(dim(SAel_all_max))
  }
}

SAel_all_max1 = as.data.frame(SAel_all_max)
dim(SAel_all_max1)
head(SAel_all_max1)
table(SAel_all_max1$y_true)

# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------
SAel_all_max1$agg_score = rowMeans(SAel_all_max1[, startsWith(colnames(SAel_all_max1), "y_pred_MA")])
length(SAel_all_max1$agg_score)
summary(SAel_all_max1$agg_score)

SAel_all_max1$pep_length = nchar(str_replace_all(SAel_all_max1$peptide, "X",""))
head(SAel_all_max1)
table(SAel_all_max1$hla, SAel_all_max1$y_true)

write.table(SAel_all_max1, paste0(dir0, "PEPPRMINT_prediction.txt"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
#system(sprintf("gzip %s", paste0(dir0, "PEPPRMINT_prediction.txt")))

# ---------------------------------------------------------------------------
# Overall Performance for PEPPRMINT
# ---------------------------------------------------------------------------
#SAel_all_max1_sub = SAel_all_max1[which(SAel_all_max1$pep_length<=11),]

#ROC
SAel_roc1 = roc.curve(scores.class0 = SAel_all_max1$agg_score, 
                 weights.class0 = SAel_all_max1$y_true, 
                 curve=TRUE)
SAel_roc1

#PPV
temp1 <- SAel_all_max1_sub[with(SAel_all_max1_sub,order(-agg_score)),]
hits = table(SAel_all_max1_sub$y_true)[2][[1]]
temp2 = temp1[1:hits,]
dim(temp2)

SAel_ppv = sum(temp2$y_true)/hits
SAel_ppv

# ---------------------------------------------------------------------------
# PEPPRMINT performance by HLA
# ---------------------------------------------------------------------------
hla = unique(SAel_all_max1$hla)

hla_auc_results = NULL 
for( h in hla){
  temp = SAel_all_max1[which(SAel_all_max1$hla == h),]
  
  #ROC
  temp_roc1 = roc.curve(scores.class0 = temp$agg_score, 
                       weights.class0 = temp$y_true, 
                       curve=TRUE)
  temp_roc1
  
  #PPV
  temp1 <- temp[with(temp,order(-agg_score)),]
  hits = table(temp$y_true)[2][[1]]
  temp2 = temp1[1:hits,]
  dim(temp2)
  
  temp_ppv = sum(temp2$y_true)/hits
  temp_ppv
  
  hla_auc_results = rbind(hla_auc_results, 
                         c(h, round(temp_roc1$auc,4), round(temp_ppv,4)))
}

hla_auc_results
hla_auc_results = as.data.frame(hla_auc_results)
colnames(hla_auc_results) = c("hla", "AUC", "PPV")
hla_auc_results$hla = as.character(hla_auc_results$hla)
hla_auc_results$AUC = as.numeric(as.character(hla_auc_results$AUC))
hla_auc_results$PPV = as.numeric(as.character(hla_auc_results$PPV))

write.table(hla_auc_results, 
            "../../../results/test_data/NetMHCpan4_1/by_hla_auc_PEPPRMINT.txt", 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


sessionInfo()
q(save="no")


