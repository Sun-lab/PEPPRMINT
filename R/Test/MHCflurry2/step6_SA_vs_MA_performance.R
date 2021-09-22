# ---------------------------------------------------------------------------
# SA vs. PEPPRMINT performance 
# ---------------------------------------------------------------------------

library(data.table)
library(PRROC)
library(stringr)

# ---------------------------------------------------------------------------
# Train data prediction
# ---------------------------------------------------------------------------
dirtrain = "../../results/PEPPRMINT/"
ma_res = "mixPep_800_bs_32_lr_0.001_e_10_layer_1_dropout_0.8_new_pi_weight_0.1_decriter_2_split0_Aug2_iter10_train_pred_all.txt.gz"
ma = fread(paste0(dirtrain, ma_res))
head(ma)
ma$pep_line = paste(ma$peptide, ma$sample, sep = ";")
ma_max_split = ma[, .(max(y_pred), max(y_pred_pMHC)), by = .(y_true, pep_line)]
dim(ma_max_split)
colnames(ma_max_split)[3:4] = c("y_pred_MA","y_pred_SA")

# MA data only 
ma_max_split$cell_line = sub(".*;", "", ma_max_split$pep_line) 
ma_max_split1 = ma_max_split[which(substr(ma_max_split$cell_line,1,3) != "HLA"),]
table(ma_max_split1$cell_line)
dim(ma_max_split1)

#Add # HLA in cell-line 
allele.list = read.delim(sprintf("%sallelelist.txt", 
                                 "../../data/NetMHCpan4_1_train/"), 
                         header=FALSE, sep = "", fill = TRUE, as.is = TRUE)
dim(allele.list)
allele.list = unique(allele.list)
dim(allele.list)

names(allele.list) = c("cell_line", "hla")
head(allele.list)

hla.list = strsplit(allele.list$hla, split=",")
n.hla    = sapply(hla.list, length)
table(n.hla)
allele.list$n.hla = n.hla

ma_max_split2 = merge(ma_max_split1, allele.list, by = "cell_line", 
                      sort = FALSE)
dim(ma_max_split2)
dim(ma_max_split1)
head(ma_max_split2)
table(ma_max_split2$n.hla)


for(i in min(ma_max_split2$n.hla):max(ma_max_split2$n.hla)){
  tempds = ma_max_split2[which(ma_max_split2$n.hla==i),]
  
  roc_ma_nhla = roc.curve(scores.class0 = tempds$y_pred_MA, 
                           weights.class0 = tempds$y_true, 
                           curve=TRUE)
  roc_ma_nhla
  
  roc_sa_nhla = roc.curve(scores.class0 = tempds$y_pred_SA, 
                           weights.class0 = tempds$y_true, 
                           curve=TRUE)
  roc_sa_nhla
  print(paste0(i, " HLA per cell line: MA AUC = ", round(roc_ma_nhla$auc,3), 
              ", SA AUC = ", round(roc_sa_nhla$auc, 3)))
}

# AUC
roc_ma_split = roc.curve(scores.class0 = ma_max_split1$y_pred_MA, 
                         weights.class0 = ma_max_split1$y_true, 
                         curve=TRUE)
roc_ma_split

roc_sa_split = roc.curve(scores.class0 = ma_max_split1$y_pred_SA, 
                         weights.class0 = ma_max_split1$y_true, 
                         curve=TRUE)
roc_sa_split

# PPV
temp1 <- ma_max_split[with(ma_max_split,order(-y_pred)),]
hits = table(ma_max_split$y_true)[2][[1]]
temp2 = temp1[1:hits,]
dim(temp2)

ppv_sa = sum(temp2$y_true)/hits
ppv_sa


# ---------------------------------------------------------------------------
# SA ROC - MHCflurry-2.0 test data prediction
# ---------------------------------------------------------------------------
sa_res = "pMHCpan_MHCflurry2_800_nocw_full_Jul26_test_pred_all.txt.gz"
sa = fread(paste0(dir0, sa_res))
sa_max = sa[, .(max(y_pred_mix)), by = .(y_true, peptide)]
dim(sa_max)
colnames(sa_max)[3] = "y_pred"

# AUC
roc_sa = roc.curve(scores.class0 = sa_max$y_pred, weights.class0 = sa_max$y_true, 
                   curve=TRUE)
roc_sa

# PPV
temp1 <- sa_max[with(sa_max,order(-y_pred)),]
hits = table(sa_max$y_true)[2][[1]]
temp2 = temp1[1:hits,]
dim(temp2)

ppv_sa = sum(temp2$y_true)/hits
ppv_sa

# ---------------------------------------------------------------------------
# aggregate 4 MA ROC MHCflurry-2.0 results 
# ---------------------------------------------------------------------------

# Extract files 
dir0 = "../../results/test_data/MHCflurry2/"
fall = list.files(dir0, 
                  pattern = "test_pred_all.txt.gz")
fall = fall[grep("800",fall)]
fall = fall[grep("MA",fall)]
fall = fall[-grep("split4",fall)]
fall
#combine datasets across splits (after taking max of each peptide)
ma_max = NULL 
for(f in fall){
  print("*-------------------------------------------------------------------------*")
  it = str_split(f, "_")[[1]][22]
  print(it)
  p_pred = fread(paste0(dir0, f))
  
  dim(p_pred)
  head(p_pred)
  
  d50 = p_pred[, .(max(y_pred_mix)), by = .(y_true, peptide)]
  colnames(d50)[3] = "y_pred_mix"
  head(d50)
  print(dim(d50))
  # Note: this dimension is correct if remove duplicates within cell line 
  
  #create cell line
  d50$cell_line = sub("\\;.*", "", d50$peptide)
  head(d50)
  
  # ---------------------------------------------------------------------------
  # Combine all split/configs to one dataset
  # ---------------------------------------------------------------------------
  colnames(d50)[which(colnames(d50)=="y_pred_mix")] = paste("y_pred_mix",
                                                            str_split(f, "_")[[1]][3],
                                                            str_split(f, "_")[[1]][22], 
                                                            sep = "_")
  colnames(d50)[which(colnames(d50)=="hla")] = paste("hla",
                                                            str_split(f, "_")[[1]][3],
                                                            str_split(f, "_")[[1]][22], 
                                                            sep = "_")
  if(str_split(f, "_")[[1]][22]=="split0"){
    ma_max = d50
    print(dim(ma_max))
  }else{
    ma_max = merge(ma_max, d50, 
                    by = c("peptide", "y_true", "cell_line"))
    print(dim(ma_max))
  }
}
head(ma_max)
table(ma_max$y_true)
ma_max = as.data.frame(ma_max)

#aggregate across 4 splits
ma_max$ag_score = rowMeans(ma_max[, startsWith(colnames(ma_max), "y_pred_mix")])

# calculate ROC curve 
roc_ma = roc.curve(scores.class0 = ma_max$ag_score, weights.class0 = ma_max$y_true, 
                 curve=TRUE)
roc_ma

# PPV

# PPV
temp1 <- ma_max[with(ma_max,order(-ag_score)),]
hits = table(ma_max$y_true)[2][[1]]
temp2 = temp1[1:hits,]
dim(temp2)

ppv_ma = sum(temp2$y_true)/hits
ppv_ma




