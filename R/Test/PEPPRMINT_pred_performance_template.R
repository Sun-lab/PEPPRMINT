# ---------------------------------------------------------------------------
# Analyze PEPPRMINT prediction results 
# ---------------------------------------------------------------------------
library(data.table)
library(PRROC)
library(stringr)

# set directory where your prediction results are located to dir0
dir0 = "../../results/test_data/MHCflurry2/"

#set the file name for aggregated score
PEPPRMINT_file_name = "PEPPRMINT_MHCflurry2_pred_Aug5.txt"
# ---------------------------------------------------------------------------
# Import data 
# ---------------------------------------------------------------------------

fall = list.files(dir0, 
                  pattern = "test_pred_all.txt.gz")
fall= fall[grep("MA", fall)]
#check that you have 15 prediction files corresponding to each model
fall

# ---------------------------------------------------------------------------
#extract cell_line (20)
# ---------------------------------------------------------------------------
p =  fread(paste0(dir0, fall[1]))
p50 = p[, .(max(y_pred_mix)), by = .(y_true, peptide)]

p50$cell_line = sub("\\;.*", "", p50$peptide)
head(p50)
cell_line = unique(p50$cell_line)
cell_line

hits_tab  = table( p50$cell_line, p50$y_true)
hits_tab

hits_tab1 = hits_tab[,2]
hits_tab1


# ---------------------------------------------------------------------------
#function for each iteration
# ---------------------------------------------------------------------------
all_max = NULL 

for(f in fall){
  print("*-------------------------------------------------------------------------*")
  it = str_split(f, "_")[[1]][3]
  #print(it)
  
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
  # Overall split/config performance (validate python results)
  # ---------------------------------------------------------------------------
  # calculate ROC curve 
  roc1 = roc.curve(scores.class0 = d50$y_pred_mix, weights.class0 = d50$y_true, 
                   curve=TRUE)
  roc1
  #plot(roc1)
  
  #pr1  = pr.curve(scores.class0 = d50$y_pred_mix, weights.class0 = d50$y_true,
  #                curve=TRUE)
  #pr1
  #plot(pr1)
  
  # PPV
  temp1 <- d50[with(d50,order(-y_pred_mix)),]
  hits = table(d50$y_true)[2][[1]]
  temp2 = temp1[1:hits,]
  dim(temp2)
  
  ppv = sum(temp2$y_true)/hits
  ppv
  
  print(paste0("AUC, PPV for config: ", f))
  print(paste("AUC =", round(roc1$auc,3)," PPV = ", round(ppv,3), sep = " "))
  
  # ---------------------------------------------------------------------------
  # Combine all split/configs to one dataset
  # ---------------------------------------------------------------------------
  colnames(d50)[which(colnames(d50)=="y_pred_mix")] = paste("y_pred_mix",
                                                            str_split(f, "_")[[1]][2],
                                                            str_split(f, "_")[[1]][3], 
                                                            sep = "_")
  if(str_split(f, "_")[[1]][2]=="200" & str_split(f, "_")[[1]][3]=="split0"){
    all_max = d50
    print(dim(all_max))
  }else{
    all_max = merge(all_max, d50, 
                    by = c("peptide", "cell_line", "y_true"))
    print(dim(all_max))
  }
}

all_max1 = as.data.frame(all_max)
dim(all_max1)
head(all_max1)
table(all_max1$y_true)

# ---------------------------------------------------------------------------
# Aggregate acroos 15 models
# ---------------------------------------------------------------------------
all_max1$cv_score = rowMeans(all_max1[, startsWith(colnames(all_max1), "y_pred_mix")])
length(all_max1$cv_score)
summary(all_max1$cv_score)

summary(all_max1$cv_score[which(all_max1$y_true ==1)])
summary(all_max1$cv_score[which(all_max1$y_true ==0)])

write.table(all_max1, paste0(dir0, PEPPRMINT_file_name),
                             col.names = TRUE, row.names = FALSE, 
                             quote = FALSE, sep = "\t")

# ---------------------------------------------------------------------------
# Overall Performance for AGGREGATE SCORES
# ---------------------------------------------------------------------------
#ROC
roc1 = roc.curve(scores.class0 = all_max1$cv_score, 
                 weights.class0 = all_max1$y_true, 
                 curve=TRUE)
roc1
plot(roc1)

#PPV
temp1 <- all_max1[with(all_max1,order(-cv_score)),]
hits = table(all_max1$y_true)[2][[1]]
temp2 = temp1[1:hits,]
dim(temp2)

ppv = sum(temp2$y_true)/hits
ppv

# ---------------------------------------------------------------------------
# AGGREGATE SCORES: Cell_line performance
# ---------------------------------------------------------------------------
cl_auc_results = NULL
cell_lines = unique(all_max1$cell_line)
length(cell_lines)

for(c in cell_lines){
  print(c)
  temp = all_max1[which(all_max1$cell_line== c), ]
  #print(paste(c, hits, sep = ": "))
  
  #AUC 
  roc1 = roc.curve(scores.class0 = temp$cv_score, weights.class0 = temp$y_true, 
                   curve=TRUE)
  roc1
  
  # PPV
  temp1 <- temp[with(temp,order(-cv_score)),]
  hits = table(temp$y_true)[2][[1]]
  temp2 = temp1[1:hits,]
  print(hits)
  dim(temp2)
  
  ppv = sum(temp2$y_true)/hits
  ppv
  
  cl_auc_results = rbind(cl_auc_results, 
                         c(c, round(roc1$auc,4), round(ppv,4)))
  
}

cl_auc_results
cl_auc_results = as.data.frame(cl_auc_results)
colnames(cl_auc_results) = c("sample_id", 'AUC', 'PPV')
cl_auc_results$sample_id = as.character(cl_auc_results$sample_id)
cl_auc_results$AUC = as.numeric(as.character(cl_auc_results$AUC))
cl_auc_results$PPV = as.numeric(as.character(cl_auc_results$PPV))


sessionInfo()
q(save="no")

