#---------------------------------------------------------------------------------
# Modify MA data so that the number of negative peptides is 
# of the same ratio to the number of positive peptides
#---------------------------------------------------------------------------------

library(data.table)
library(stringr)

dir0 = "../../../data/NetMHCpan4_1_train"
dir1 = "../../../data/MA_data"

for(i in 0:4){
  fnm   = sprintf("%s/train_v4_el_multi_HLA_%d.txt.gz", dir1, i)
  cat(i, " ", fnm, "\n")
  
  dat_i = fread(fnm)
  dim(dat_i)
  dat_i[1:5,]
  
  t2 = tapply(dat_i$peptide, dat_i$cell_line, anyDuplicated)
  stopifnot(all(t2 == 0))
  
  stopifnot(all(str_length(dat_i$peptide) == 15))
  dat_i$len = 15 - str_count(dat_i$peptide, "X")
  
  tb1 = table(dat_i$len)
  tb2 = table(dat_i$binder, dat_i$len)
  
  print(tb1)
  print(tb2)
  print(tb2[1,]/tb2[2,])
  
  # keep the number for 9aa peptide and select 5x negatives 
  # for all other lengths
  
  w2kp = NULL
  
  set.seed(111)
  
  for(l1 in c(8, 10:15)){
    print(l1)
    w_pos  = which(dat_i$len == l1 & dat_i$binder == 1)
    w_neg  = which(dat_i$len == l1 & dat_i$binder == 0)
    n_pos  = length(w_pos)
    n_neg  = 5*n_pos
    
    w2kp   = c(w2kp, w_pos, sample(w_neg, n_neg))
  }
  
  dat_i_new = dat_i[w2kp,]
  dat_i_new = rbind(dat_i_new, dat_i[which(dat_i$len ==9),])
  dim(dat_i_new)
  
  tb1 = table(dat_i_new$len)
  tb2 = table(dat_i_new$binder, dat_i_new$len)
  
  print(tb1)
  print(tb2)
  print(tb2[1,]/tb2[2,])
  
  fnm1_new = sub("multi_HLA_", "multi_HLA_balanced_", fnm)
  fnm1_new = sub(".gz", "", fnm1_new)
  
  cols2kp = c("peptide", "binder", "cell_line")
  
  write.table(dat_i_new[,..cols2kp], fnm1_new, sep="\t", quote = FALSE, 
              row.names = FALSE)
  system(sprintf("gzip %s", fnm1_new))
}

sessionInfo()
q(save="no")
