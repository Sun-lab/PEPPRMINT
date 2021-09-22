# ---------------------------------------------------------------------------
# Analyze PEPPRMINT prediction results 
# ---------------------------------------------------------------------------
library(data.table)
library(PRROC)
library(stringr)

# ---------------------------------------------------------------------------
# Import data 
# ---------------------------------------------------------------------------
dir0 = "../../../results/test_data/MHCflurry2/"
fall = list.files(dir0, 
                  pattern = "test_pred_all.txt.gz")
fall= fall[grep("MA", fall)]
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
  it = str_split(f, "_")[[1]][23]
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
                                                            str_split(f, "_")[[1]][3],
                                                            str_split(f, "_")[[1]][22], 
                                                            sep = "_")
  if(str_split(f, "_")[[1]][3]=="200" & str_split(f, "_")[[1]][22]=="split0"){
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
# Aggregate
# ---------------------------------------------------------------------------
all_max1$cv_score = rowMeans(all_max1[, startsWith(colnames(all_max1), "y_pred_mix")])
length(all_max1$cv_score)
summary(all_max1$cv_score)

summary(all_max1$cv_score[which(all_max1$y_true ==1)])
summary(all_max1$cv_score[which(all_max1$y_true ==0)])

write.table(all_max1, paste0(dir0,"PEPPRMINT_MHCflurry2_pred_Aug5.txt"),
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

write.table(cl_auc_results, 
            "../../../results/test_data/MHCflurry2/cell_line_performance_PEPPRMINT.txt", 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)



sessionInfo()
q(save="no")




# ---------------------------------------------------------------------------
#import perforance results 
# ---------------------------------------------------------------------------
comp= read.delim("../../../../data/test_data/MHCflurry2/performance_results.txt", 
                 sep = "\t")


comp 

iter10_res = all_iter[,5:8]
iter10_res

iter10_res$sample_id = str_replace(iter10_res$sample, "_iter10", "")
head(iter10_res)


comp_fin = merge(comp, iter10_res, by = "sample_id")
comp_fin$betterNetEL = ifelse(comp_fin$`mix_Aug2 AUC` > comp_fin$AUC..NetMHCpan.4.0.EL, "yes", "no")
comp_fin$betterNetBA = ifelse(comp_fin$`mix_Aug2 AUC` > comp_fin$AUC..NetMHCpan.4.0.BA, "yes", "no")
comp_fin$betterMHCf2 = ifelse(comp_fin$`mix_Aug2 AUC` > comp_fin$AUC..MHCflurry.2.0.BA, "yes", "no")
table(comp_fin$betterNetEL)
table(comp_fin$betterNetBA)
table(comp_fin$betterMHCf2)

comp_fin$betterNetELppv = ifelse(comp_fin$`mix_Aug2 PPV` > comp_fin$PPV..NetMHCpan.4.0.EL, "yes", "no")
comp_fin$betterNetBAppv = ifelse(comp_fin$`mix_Aug2 PPV` > comp_fin$PPV..NetMHCpan.4.0.BA, "yes", "no")
comp_fin$betterMHCf2ppv = ifelse(comp_fin$`mix_Aug2 PPV` > comp_fin$PPV..MHCflurry.2.0.BA, "yes", "no")
table(comp_fin$betterNetELppv)
table(comp_fin$betterNetBAppv)
table(comp_fin$betterMHCf2ppv)

table(comp_fin$betterNetELppv, comp_fin$bette)

write.table(comp_fin[c("sample_id","HITS", "AUC..NetMHCpan.4.0.EL","AUC..NetMHCpan.4.0.BA","AUC..MHCflurry.2.0.BA",   
           "mix_Aug2 AUC","AUC..MixMHCpred.2.0.2","PPV..NetMHCpan.4.0.EL","PPV..NetMHCpan.4.0.BA",
           "PPV..MHCflurry.2.0.BA", "mix_Aug2 PPV", "PPV..MixMHCpred.2.0.2")], 
           "../../../../results/test_data/MHCflurry2/comparison_iteration_mixpep800_Aug2.txt", 
           sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)






# Try for one dataset 
p_pred = fread("../../A_Dissertation/A_GitHubDesktop/mixPep/results/test_data/")

dim(p_pred)
head(p_pred)


d50 = p_pred[, .(max(y_pred_mix)), by = .(y_true, peptide)]
colnames(d50)[3] = "y_pred"
head(d50)
dim(d50)
# Note: this dimension is correct if remove duplicates within cell line 

# ---------------------------------------------------------------------------
# Overall performance (validate python results)
# ---------------------------------------------------------------------------
# calculate ROC curve 
roc1 = roc.curve(scores.class0 = d50$y_pred, weights.class0 = d50$y_true, 
                 curve=TRUE)
roc1
plot(roc1)

pr1  = pr.curve(scores.class0 = d50$y_pred, weights.class0 = d50$y_true,
                curve=TRUE)
pr1
plot(pr1)
# ---------------------------------------------------------------------------
#extract cell_line
# ---------------------------------------------------------------------------
d50$cell_line = sub("\\;.*", "", d50$peptide)
head(d50)
cell_line = unique(d50$cell_line)
cell_line

# ---------------------------------------------------------------------------
# Performance for each Cell Line (compare with MHCflurry-2.0 comparison results)
# ---------------------------------------------------------------------------
p_auc_results = NULL
for(c in cell_line){
  temp = d50[which(d50$cell_line== c), ]
  
  roc1 = roc.curve(scores.class0 = temp$y_pred, weights.class0 = temp$y_true, 
                   curve=TRUE)
  roc1
  
  
  pr1  = pr.curve(scores.class0 = temp$y_pred, weights.class0 = temp$y_true,
                  curve=TRUE)
  pr1
  
  p_auc_results = rbind(p_auc_results, c(c, round(roc1$auc,4), 
                                         round(pr1$auc.integral,4)))
  
  #plot(roc1, main = paste("ROC", c, sep = " "))
  #plot(pr1, main = paste("PR", c, sep = " "))
  
}

p_auc_results = as.data.frame(p_auc_results)
colnames(p_auc_results) = c("sample_id", 'pM_Jul26 ROC', "pM_Jul26 PR" )
p_auc_results$`pM_Jul26 ROC` = as.numeric(as.character(p_auc_results$`pM_Jul26 ROC`))
p_auc_results$`pM_Jul26 PR` = as.numeric(as.character(p_auc_results$`pM_Jul26 PR`))
p_auc_results


# ---------------------------------------------------------------------------
# Subset to peptides not in NetMHCpan-4.1 training data
# ---------------------------------------------------------------------------
nodup = fread("../../../data/test_data/MHCflurry2/MHCflurry2_test_el_multi_HLA_no_overlap.txt", 
              sep = "\t")
head(nodup)
nodup = nodup[,c("peptide", "y_true")]

all_max1_bind = all_max1[which(all_max1$y_true ==1),]
all_max1_bind_nodup = merge(all_max1_bind, nodup, by = intersect(colnames(all_max1_bind), 
                                                                 colnames(nodup)), 
                            sort = FALSE)
dim(nodup)
dim(all_max1_bind_nodup)

all_max1_noover = rbind(all_max1_bind_nodup, all_max1[which(all_max1$y_true ==0),])
dim(all_max1_noover)

#ROC
roc1noover = roc.curve(scores.class0 = all_max1_noover$cv_score, 
                       weights.class0 = all_max1_noover$y_true, 
                       curve=TRUE)
roc1noover
plot(roc1oover)

#PPV
temp1 <- all_max1_noover[with(all_max1_noover,order(-cv_score)),]
hits = table(all_max1_noover$y_true)[2][[1]]
temp2 = temp1[1:hits,]
dim(temp2)

ppvnoover = sum(temp2$y_true)/hits
ppvnoover

