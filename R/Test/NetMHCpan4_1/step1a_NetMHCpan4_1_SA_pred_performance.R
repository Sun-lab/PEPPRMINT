# ---------------------------------------------------------------------------
# Performance of NetMHCpan-4.1 on NetMHCpan-4.1 SAel Testing Data
# ---------------------------------------------------------------------------
library(ggplot2)
library(pROC)
library(PRROC)
library(data.table)
library(stringr)

# ---------------------------------------------------------------------------
# Import predictins
# ---------------------------------------------------------------------------
#fall = list.files("../../../../netMHCpan-4.1/test/_predictions/", pattern=".xls")
fall = list.files("../../../results/test_data/NetMHCpan4_1/_predictions/", pattern=".xls")
fall
dir0 = "../../../data/test_data/NetMHCpan4_1/"

SAel_net_pred = NULL 
for(file in fall){
  print("--------------------")
  print(file)
  #import data set
  ds = fread(paste("../../../results/test_data/NetMHCpan4_1/_predictions/", 
                   file, sep = ""), skip = 1, 
             header = TRUE)
  head(ds)
  print(dim(ds))
  
  
  #add HLA column to dataset
  ds$HLA = str_extract(file, "[^_]+")
  
  # Remove duplicated peptides in HLA
  ds$pep_line = paste(ds$HLA, ds$Peptide, sep = ";")
  ds1= ds[!duplicated(ds$pep_line),]
  
  #add true binding status
  file_true = str_replace(str_extract(file, "[^_]+"), ":", "_")
  truds = fread(paste(dir0, 
                     file_true, ".txt", sep = ""), header = FALSE)
  colnames(truds) = c("Peptide", "binder", "HLA")
  head(truds)
  print("Binding status dataset ")
  print(dim(truds))
  truds$pep_line = paste(truds$HLA, truds$Peptide, sep = ";")
    # remove duplicated lines
  truds1= truds[!duplicated(truds$pep_line),]
  print("Binding status dataset after removing duplciated lines")
  print(dim(truds1))
  
  
  # Merge datasets: add true bind to prediction
  ds1 = merge(ds1, truds1, by = intersect(colnames(ds1), colnames(truds1)), 
              sort = FALSE)
  print("merged dataset")
  print(dim(ds1))
  head(ds1)
  
  SAel_net_pred = rbind(SAel_net_pred, ds1)
  
  print(dim(SAel_net_pred))
  
}

head(SAel_net_pred)
table(SAel_net_pred$HLA)
dim(SAel_net_pred)
table(SAel_net_pred$binder)

write.table(SAel_net_pred, 
            "../../../results/test_data/NetMHCpan4_1/NetMHCpan4_1_predictions.txt", 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# ---------------------------------------------------------------------------
# AUC analysis 
# ---------------------------------------------------------------------------
rocnet_SAel = roc.curve(scores.class0 = SAel_net_pred$`EL-score`, 
                weights.class0 = SAel_net_pred$binder)

rocnet_SAel

#PPV
temp1 <- SAel_net_pred[with(SAel_net_pred,order(-`EL-score`)),]
hits = table(SAel_net_pred$binder)[2][[1]]
temp2 = temp1[1:hits,]
dim(temp2)

SAel_ppv = sum(temp2$binder)/hits
SAel_ppv

# ---------------------------------------------------------------------------
# NetMHCpan-4.1 performance by HLA
# ---------------------------------------------------------------------------
hla = unique(SAel_net_pred$HLA)

hla_auc_results_net = NULL 
for( h in hla){
  print(h)
  temp = SAel_net_pred[which(SAel_net_pred$HLA == h),]
  
  #ROC
  temp_roc1 = roc.curve(scores.class0 = temp$`EL-score`, 
                        weights.class0 = temp$binder, 
                        curve=TRUE)
  temp_roc1
  
  #PPV
  temp1 <- temp[with(temp,order(-`EL-score`)),]
  hits = table(temp$binder)[2][[1]]
  temp2 = temp1[1:hits,]
  dim(temp2)
  
  temp_ppv = sum(temp2$binder)/hits
  temp_ppv
  
  hla_auc_results_net = rbind(hla_auc_results_net, 
                          c(h, round(temp_roc1$auc,4), round(temp_ppv,4)))
}

hla_auc_results_net
hla_auc_results_net = as.data.frame(hla_auc_results_net)
colnames(hla_auc_results_net) = c("hla", "AUC_net", "PPV_net")
hla_auc_results_net$hla = as.character(hla_auc_results$hla)
hla_auc_results_net$AUC_net = as.numeric(as.character(hla_auc_results_net$AUC_net))
hla_auc_results_net$PPV_net = as.numeric(as.character(hla_auc_results_net$PPV_net))

write.table(hla_auc_results_net, 
            "../../../results/test_data/NetMHCpan4_1/by_hla_auc_NetMHCpan4_1.txt", 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

sessionInfo()
q(save="no")

# ---------------------------------------------------------------------------
# Density Plots
# ---------------------------------------------------------------------------
net_pred$true_bind = ifelse(net_pred$binder==1, "binder", "nonbinder")
p3 = ggplot(net_pred, aes(x = `EL-score`, fill = true_bind)) + 
  geom_density(alpha = 0.5) + 
  ggtitle(paste("Density plot for test NetMHCpan-4.1 predictions 
                all peptides: AUC = ", round(roc$auc,3), ", PR= ", round(pr$auc.integral,3), 
                sep = "")) + 
  xlab("NetMHCpan-4.1 prediction")

p3_1 = ggplot(net_pred[which(net_pred$binder==1),], 
              aes(x = `EL-score`, fill = true_bind)) + 
  geom_density(alpha = 0.5) + 
  ggtitle(paste("Density plot for test NetMHCpan-4.1 predictions 
                Binders: AUC = ", round(roc$auc,3), ", PR= ", round(pr$auc.integral,3), 
                sep = "")) + 
  xlab("NetMHCpan-4.1 prediction")

p3_2 = ggplot(net_pred[which(net_pred$binder==0),], 
              aes(x = `EL-score`, fill = true_bind)) + 
  geom_density(alpha = 0.5) + 
  ggtitle(paste("Density plot for test NetMHCpan-4.1 predictions 
                Nonbinders: AUC = ", round(roc$auc,3), ", PR= ", round(pr$auc.integral,3), 
                sep = "")) + 
  xlab("NetMHCpan-4.1 prediction")


pdf(paste("../../../../figures/mixPep/NetMHCpan_4_1_test_prediction_Jun7.pdf", sep = ""), 
    width=8.5, height=5,onefile=T)

print(p3)
print(p3_1)
print(p3_2)

dev.off()

