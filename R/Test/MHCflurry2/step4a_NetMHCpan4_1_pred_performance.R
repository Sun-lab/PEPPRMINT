# ------------------------------------------------------------------------------------
# NetMHCpan-4.1 predictions in MHCflurry-2.0
# ------------------------------------------------------------------------------------
library(stringr)
library(data.table)
library(matrixStats)
library(dplyr)
library(PRROC)
dir0 = #set to the directory where your NetMHCpan-4.1 predictions are 
#"../../../../netMHCpan-4.1/test/_predictions/"
# ------------------------------------------------------------------------------------
# Import true binding status
# ------------------------------------------------------------------------------------
mhcflurry= fread("../../../data/test_data/MHCflurry2/MHCflurry2_test_el_multi_HLA.txt.gz")
unique(mhcflurry$cell_line)
mhcflurry$cell_line_rpl = str_replace_all( mhcflurry$cell_line, "/","_")
head(mhcflurry)
dim(mhcflurry)

# ------------------------------------------------------------------------------------
# Find Prediction Files 
# ------------------------------------------------------------------------------------
fall = list.files(dir0, pattern="MHCflurry2_Net")
fall = fall[grep("xls", fall)]
fall

# ------------------------------------------------------------------------------------
# Caculate AUC and PPV per cell line and combine all cell-lines
# ------------------------------------------------------------------------------------
net_pred = NULL 
cl_net_auc = NULL 
for(f in fall){
  print(f)
  tempds = read.csv(paste(dir0, f, sep = ""), header = TRUE, 
                    sep = "\t", skip = 1)
  print(nrow(tempds))
  
  tempds1 = select(tempds,contains("EL.score"))
  tempds1 = as.matrix(tempds1)
  
  tempds$max_pred = rowMaxs(tempds1)
 
  #add cell_line to dataset
  c = str_match(f, "MHCflurry2_Net_\\s*(.*?)\\s*_pred.xls")[,2]
  tempds$cell_line <- c
  
  tempds_max = tempds[,c("Pos", "Peptide", "ID", "core", "icore", 
                         "Ave", "NB", "max_pred", "cell_line")]
  head(tempds_max)
  
  # Remove peptides that are duplicated within the same cell-line
  tempds_max$pep_line = paste(tempds_max$cell_line, tempds_max$Peptide, sep = ";")
  tempds_max1 = tempds_max[!duplicated(tempds_max$pep_line),]
  dim(tempds_max1)
  
  #add true binding status to dataset
  mhctemp = mhcflurry[which(mhcflurry$cell_line_rpl==c),]
  mhctemp$pep_line = paste(mhctemp$cell_line_rpl, 
                           str_replace_all(mhctemp$peptide,"X", ""), 
                           sep = ";")
  mhctemp1 = mhctemp[!duplicated(mhctemp$pep_line),]
  dim(mhctemp1)
  head(mhctemp1)
  
  dim(tempds_max1)
  
  tempds_max2 = merge(tempds_max1, mhctemp1[,c("pep_line", "binder")],
                      by = "pep_line", sort = FALSE)
  head(tempds_max2)
  
  #AUC 
  roc0 = roc.curve(scores.class0 = tempds_max2$max_pred, 
                  weights.class0 = tempds_max2$binder, curve = TRUE)
  #pr0 = pr.curve(scores.class0 = tempds_max1$max_pred, 
  #              weights.class0 = tempds_max1$binder, curve = TRUE)
  
  # PPV
  hits = table(tempds_max2$binder)[2][[1]]
  temp1 <- tempds_max2[with(tempds_max2,order(-max_pred)),]
  temp2 = temp1[1:hits,]
  dim(temp2)
  
  ppv = sum(temp2$binder)/hits
  ppv
  
  cl_net_auc = rbind(cl_net_auc, 
                         c(c, round(roc0$auc,4), round(ppv,4)))
  
  print(paste0("cell_line ", c, ": ROC = ", round(roc0$auc,3), 
              ", PPV = ", round(ppv, 3)))
  
  if(f == "MHCflurry2_Net_10-002-S1-TISSUE_pred.xls"){
    net_pred = tempds_max2
  }else{
    net_pred = rbind(net_pred, tempds_max2)
  }
  print(dim(net_pred))
}

head(net_pred)
colnames(cl_net_auc) = c("sample_id", "NetMHCpan4.1_AUC", 'NetMHCpan4.1_PPV')
cl_net_auc = as.data.frame(cl_net_auc)
cl_net_auc$sample_id = as.character(cl_net_auc$sample_id)
cl_net_auc$NetMHCpan4.1_AUC = as.numeric(as.character(cl_net_auc$NetMHCpan4.1_AUC))
cl_net_auc$NetMHCpan4.1_PPV = as.numeric(as.character(cl_net_auc$NetMHCpan4.1_PPV))


write.table(net_pred,
            "../../results/test_data/MHCflurry2/NetMHCpan4_1_MHCflurry2_test_pred.txt", 
            row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = "\t")

write.table(cl_net_auc,
            "../../results/test_data/MHCflurry2/NetMHCpan4_1_MHCflurry2_cell_line_performance.txt", 
            row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = "\t")

# ------------------------------------------------------------------------------------
# OVERALL AUC and PPV
# ------------------------------------------------------------------------------------
netroc = roc.curve(scores.class0 = net_pred$max_pred, 
                   weights.class0 = net_pred$binder, curve = TRUE)

netroc

# PPV
hits = table(net_pred$binder)[2][[1]]
temp1 <- net_pred[with(net_pred,order(-max_pred)),]
temp2 = temp1[1:hits,]
dim(temp2)

ppvnet = sum(temp2$binder)/hits
ppvnet


sessionInfo()
q(save="no")


