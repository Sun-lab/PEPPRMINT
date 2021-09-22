#---------------------------------------------------------------------------------
# Step 1. Extract SA data for Training pMHCpan
# Note: For SA HLA-I EL peptide data from NetMHCpan-4.1
#   Steps: 
#       1. Check data from NetMHCpan-4.1 (NetMHCpan_train)
#       2. Extract EL data that are associated with a single HLA allele
#---------------------------------------------------------------------------------

library(data.table)
library(stringr)

dir0 = "../../data/NetMHCpan4_1_train/"

#---------------------------------------------------------------------------------
# Check that SA EL training data is the same as in NetMHCpan-4.1
# # HLA = 142
# # binders = 218,962
# # nonbinders = 3,813,877
# Note: for NetMHCpan-4.1 data, the validation set is not provided. 
#---------------------------------------------------------------------------------
# import 5 splits of training data 

c0 = fread(paste0(dir0, "c000_el"))
c1 = fread(paste0(dir0, "c001_el"))
c2 = fread(paste0(dir0, "c002_el"))
c3 = fread(paste0(dir0, "c003_el"))
c4 = fread(paste0(dir0, "c004_el"))

dim(c0)
c0[1:2,]

# Combine all splits of data 
dat0 = rbind(c0, c1)
dat0 = rbind(dat0, c2)
dat0 = rbind(dat0, c3)
dat0 = rbind(dat0, c4)
dim(dat0)

head(dat0)
colnames(dat0) = c("peptide", "binder", "HLA")

table(dat0$binder)
sort(unique(dat0$HLA))

# --------------------------------------------------------------------------
# Extract EL data that are associated with a single HLA allele
# two ways: 
    # 1: using allelelist and choosign lines with only 1 HLA
    # 2: using Supplementary Table 4
# --------------------------------------------------------------------------
# using Allele
# there are some duplicated rows
allele.list = read.delim(sprintf("%sallelelist.txt", dir0), header=FALSE, 
                         sep = "", fill = TRUE, as.is = TRUE)
dim(allele.list)
allele.list = unique(allele.list)
dim(allele.list)

names(allele.list) = c("cell_line", "hla")
head(allele.list)

allele.list[which(allele.list$cell_line == "Line.1"),]

hla.list = strsplit(allele.list$hla, split=",")
n.hla    = sapply(hla.list, length)
table(n.hla)

alleles = allele.list$cell_line[which(n.hla ==1)]
head(alleles)
length(alleles)

# Matching SA lines with Supplementary Table 4.
# 8 H-2, 3 DLA, 1 Mamu, A, B, C
alleles0 = unique(dat0$HLA)
length(alleles0)
head(alleles0)

table(alleles0 %in% allele.list$cell_line)
table(alleles0 %in% alleles)

alleles1 = alleles0[grep("HLA-", alleles0, fixed=TRUE)]
length(alleles1)

setdiff(intersect(alleles0, alleles), alleles1)

alleles_sup3 = alleles1
table(alleles_sup3 %in% allele.list$cell_line)

alleles_sa   = allele.list[which(allele.list$cell_line %in% alleles_sup3),]
dim(alleles_sa)
head(alleles_sa)
table(alleles_sa$cell_line == alleles_sa$hla)

write.table(alleles_sa, "../../data/SA_data/alleles_SAonly.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

length(alleles_sup3)
sort(table(dat0$HLA[which(dat0$HLA %in% alleles_sup3)]))

dat1 = dat0[which(dat0$HLA %in% alleles_sup3),]
length(unique(dat1$HLA))
dim(dat1)
head(dat1)

table(dat1$binder)
tb1 = table(dat1$binder, dat1$HLA)
dim(tb1)
tb1[,1:5]

print("ratio of non-binders to binders")
summary(tb1[1,]/tb1[2,])

dat1$pep_length = nchar(dat1$peptide)
table(dat1$pep_length, dat1$binder)
rm(dat1)
gc()

set.seed(123)

# --------------------------------------------------------------------------
# format and save training data for each split
# --------------------------------------------------------------------------

for(i in 0:4){
  set.seed(101)
  
  train.data0 = fread(sprintf("%sc00%d_el", dir0, i))

  print("-----------------------------------------------------")
  print(sprintf("i=%d", i))
  print("-----------------------------------------------------")
  
  names(train.data0)  = c("peptide", "binder", "HLA")
  
  print(table(alleles_sup3 %in% train.data0$HLA))
  
  # subset data to SA alleles only 
  train.data0 = train.data0[which(train.data0$HLA %in% alleles_sup3),]
  table(train.data0$HLA)
  table(train.data0$binder)

  t1 = table(str_length(train.data0$peptide), train.data0$binder)
  t1 = cbind(t1, t1[,2]/rowSums(t1))
  t1
  
  ww1 = train.data0$HLA %in% alleles_sup3 & 
    str_length(train.data0$peptide)<=11 &
    str_length(train.data0$peptide)>=8

  print(paste("removed", toString(sum(str_length(train.data0$peptide)<8)), 
              "and", toString(sum(str_length(train.data0$peptide)>11)), 
              "training peptides with length < 8 or > 15, respectively"))

  train.data = train.data0[which(ww1),]
  print(table(train.data$binder))
  print(table(alleles_sup3 %in% train.data$HLA))

  # check for duplication for peptide+HLA, remove all instances except original
  print("dimension of training data")
  print(dim(train.data))
  
  if(sum(duplicated(train.data[, c("peptide", "HLA")]))>0){
    temp1 = train.data[, c("peptide", "HLA")]
    train.data = train.data[!duplicated(temp1),]
    print(paste("Removed", sum(duplicated(temp1)),
                "duplicated peptide+HLA paired observations in training."))
  }else{
    print("Note: No duplicated peptide+HLA paired observations in training")
  }
  
  # create 9AA representation for training data 
  train.data$length = nchar(train.data$peptide)
  print("final data distributions by length and binder")
  print(table(train.data$length))
  print(table(train.data$binder))
  print("ratio of nonbinders:binders")
  print(round(table(train.data$binder)[1]/table(train.data$binder)[2],2))
  
  ww8  = which(train.data$length==8)
  ww9  = which(train.data$length==9)

  data_all = NULL
  data_new = NULL
  
  data8 = train.data[ww8,]
  start_pos = 1
  data8$pep_core  = paste0("X", data8$peptide)
  data8$start_pos = start_pos
  data_new = rbind(data_new, data8)
  
  for(start_pos in 2:9){
    data8 = train.data[ww8,]
    data8$pep_core  = gsub(paste0("^(.{", start_pos-1, "})(.*)$"), 
                           "\\1X\\2", data8$peptide)
    data8$start_pos = start_pos
    data_new = rbind(data_new, data8)
  }
  
  # for each non-binder, randomly choose one start position
  data_n1 = data_new[binder==1,]
  data_n0 = data_new[binder==0,]
  
  data_n0$key = paste(data_n0$peptide, data_n0$HLA, sep=":")
  stopifnot(all(table(data_n0$key) == 9))
  
  dim(data_n0)
  data_n0
  
  data_n0 = data_n0[, .SD[sample(x = .N, size = 1)], by = key]
  dim(data_n0)
  data_n0
  
  table(data_n0$start_pos)
  
  stopifnot(all(table(data_n0$key) == 1))
  data_n0[, key:=NULL]
  
  data_n1$weight = 1/9
  data_n0$weight = 1
  
  data_all = rbind(data_all, data_n1)
  data_all = rbind(data_all, data_n0)
  
  data9 = train.data[ww9,]
  data9$pep_core = data9$peptide
  data9$start_pos = 1
  data9$weight = 1
  data_all = rbind(data_all, data9)
  
  for(len1 in 10:11){
    data_l = train.data[which(train.data$length==len1),]
    len2rm = len1 - 9
    
    data_new = NULL
    
    for(k in 1:(len1 - len2rm + 1)){
      data_l$pep_core = data_l$peptide
      str_sub(data_l$pep_core, start=k, end=k+len2rm-1) = ""
      data_l$start_pos = k
      data_new = rbind(data_new, data_l)
    }
    
    data_n1 = data_new[binder==1,]
    data_n0 = data_new[binder==0,]
    
    data_n0$key = paste(data_n0$peptide, data_n0$HLA, sep=":")
    stopifnot(all(table(data_n0$key) == 10))
    
    dim(data_n0)
    data_n0
    
    data_n0 = data_n0[, .SD[sample(x = .N, size = 1)], by = key]
    dim(data_n0)
    data_n0
    table(data_n0$start_pos)
    
    stopifnot(all(table(data_n0$key) == 1))
    data_n0[, key:=NULL]
    
    data_n1$weight = 1/10
    data_n0$weight = 1
    
    data_all = rbind(data_all, data_n1)
    data_all = rbind(data_all, data_n0)
  }
  
  print("dimension of new data")
  print(dim(data_all))
  print(head(data_all))
  
  print("length of pepcore")
  print(table(nchar(data_all$pep_core)))
  
  print("table of weight vs. peptide length")
  print(table(data_all$length, data_all$weight))
  
  print("table of start_pos vs. peptide length")
  print(table(data_all$length, data_all$start_pos))
  
  # Output data files 
  dir1= "../../data/SA_data/by_split"

  fnm = sprintf("%s/train_v4_el_single_HLA_9AA_%d.txt", dir1, i)
  fwrite(data_all, fnm, sep="\t")
  system(sprintf("gzip %s", fnm))
}

gc()
sessionInfo()
q(save="no")

