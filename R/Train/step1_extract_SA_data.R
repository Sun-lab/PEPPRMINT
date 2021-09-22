#---------------------------------------------------------------------------------
# Step 1. Extract SA data for Training pMHCpan
# Note: For SA HLA-I EL peptide data from NetMHCpan-4.1
#   Steps: 
#       1. Check data from NetMHCpan-4.1 (NetMHCpan_train)
#       2. Extract EL data that are associated with a single HLA allele
#---------------------------------------------------------------------------------

library(data.table)
library(stringi)
library(ggplot2)
theme_set(theme_classic())

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
    # 1: using allelelist and choosing lines with only 1 HLA
    # 2: using Supplementary Table 4
    # 3: using human HLA only 
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

table(alleles0 %in% alleles)

alleles1 = alleles0[grep("HLA-", alleles0, fixed=TRUE)]
length(alleles1)

setdiff(intersect(alleles0, alleles), alleles1)

#alleles2 = alleles0[grep("H-2", alleles0, fixed=TRUE)]
#alleles2
#alleles3 = alleles0[grep("DLA", alleles0, fixed=TRUE)]
#alleles3
#alleles4 = alleles0[grep("Mamu", alleles0, fixed=TRUE)]
#alleles4

#alleles_sup3 = c(alleles1, alleles2, alleles3, alleles4)
alleles_sup3 = alleles1
table(alleles_sup3 %in% allele.list$cell_line)

alleles_sa   = allele.list[which(allele.list$cell_line %in% alleles_sup3),]
dim(alleles_sa)

write.table(alleles_sa, "../../data/SA_data/alleles_SAonly.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

length(alleles_sup3)
sort(table(dat0$HLA[which(dat0$HLA %in% alleles_sup3)]))

dat1 = dat0[which(dat0$HLA %in% alleles_sup3),]
length(unique(dat1$HLA))
dim(dat1)
head(dat1)

length(unique(dat1$HLA))
table(dat1$binder)
tb1 = table(dat1$binder, dat1$HLA)
dim(tb1)
tb1[,1:5]

summary(tb1[2,]/tb1[1,])

df1 = data.frame(n_sample = as.numeric(table(dat1$HLA)))
g1 = ggplot(df1, aes(x=log10(n_sample))) + xlab("log10(# of peptides per HLA)") + 
  geom_histogram(color="darkblue", fill="lightblue", bins=20)
pdf("figures/hist_n_sample_SA.pdf", width=3, height=2)
g1
dev.off()

summary(df1$n_sample)

# --------------------------------------------------------------------------
# remove duplicated peptides within the same HLA
# --------------------------------------------------------------------------

# check for duplication for peptide+HLA, remove all instances except original
if(sum(duplicated(dat1[, c("peptide", "HLA")]))>0){
  temp1 = dat1[, c("peptide", "HLA")]
  dat1  = dat1[!duplicated(temp1),]
  print(paste("Removed", sum(duplicated(temp1)),
              "duplicated peptide+HLA pairs in training.", sep = " "))
}else{
  print("Note: no duplicated peptide+HLA pairs in training" )
}

# check for peptides shared across HLA alleles

wdup = duplicated(dat1$peptide) | duplicated(dat1$peptide, fromLast=TRUE)
ndup = sum(wdup)
ndup

table(dat1[wdup,]$binder)

dat1$length = nchar(dat1$peptide)
dim(dat1)
head(dat1)

tapply(dat1$binder, dat1$length, table)

rm(dat1)

gc()

set.seed(123)

# --------------------------------------------------------------------------
# format and save training data for each split
# --------------------------------------------------------------------------
for(i in 0:4){
  train.data0 = fread(sprintf("%sc00%d_el", dir0, i))
  # test.data0  = fread(sprintf("%sc00%d_el", dir0, i))
  
  print("\n-----------------------------------------------------\n")
  print(sprintf("i=%d", i))
  print("\n-----------------------------------------------------\n")
  
  names(train.data0)  = c("peptide", "binder", "HLA")
  
  print(table(alleles_sup3 %in% train.data0$HLA))
  
  #subset data to SA alleles only 
  train.data0 = train.data0[which(train.data0$HLA %in% alleles_sup3),]
  table(train.data0$HLA)
  table(train.data0$binder)
  #print(table(alleles %in% test.data0$HLA))
  
  t1 = table(stri_length(train.data0$peptide), train.data0$binder)
  t1 = cbind(t1, t1[,2]/rowSums(t1))
  #print(t1)
  
  ww1 = train.data0$HLA %in% alleles_sup3 & 
    stri_length(train.data0$peptide)<=15 &
    stri_length(train.data0$peptide)>=8

  print(paste("removed", toString(sum(stri_length(train.data0$peptide)<8)), 
              "and", toString(sum(stri_length(train.data0$peptide)>15)), 
              "training peptides with length <8 or >15, respectively"))

  train.data = train.data0[which(ww1),]
  print(table(train.data$binder))
  print(table(alleles_sup3 %in% train.data$HLA))

  #check for duplication for peptide+HLA, remove all instances except original
  print(dim(train.data))
  if(sum(duplicated(train.data[, c("peptide", "HLA")]))>0){
    temp1 = train.data[, c("peptide", "HLA")]
    train.data = train.data[!duplicated(temp1),]
    print(paste("Removed", sum(duplicated(temp1)),
                "duplicated peptide+HLA paired observations in training."))
  }else{
    print("Note: No duplicated peptide+HLA paired observations in training")
  }
  
  #create 15 length representation for training data 
  train.data$length = nchar(train.data$peptide)
  
  print("--> Post-filtering SA distributions by length and binder <--")
  print(table(train.data$length))
  print(table(train.data$binder))
  print("--> Post-filtering SA ratio of nonbinders:binders <--")
  print(round(table(train.data$binder)[1]/table(train.data$binder)[2],2))
  
  train.data$peprep = NA
  ww8  = which(train.data$length==8)
  ww9  = which(train.data$length==9)
  ww10 = which(train.data$length==10)
  ww11 = which(train.data$length==11)
  ww12 = which(train.data$length==12)
  ww13 = which(train.data$length==13)
  ww14 = which(train.data$length==14)
  ww15 = which(train.data$length==15)
  
  train.data$peprep[ww8] = paste0(
    substr(train.data$peptide[ww8], 1, 4), "XXXXXXX", 
    substr(train.data$peptide[ww8], 5, 8))
  
  train.data$peprep[ww9] = paste0(
    substr(train.data$peptide[ww9], 1, 4), "XXX", 
    substr(train.data$peptide[ww9], 5, 5), "XXX", 
    substr(train.data$peptide[ww9], 6, 9))
  
  train.data$peprep[ww10] = paste0(
    substr(train.data$peptide[ww10], 1, 4), "XXX", 
    substr(train.data$peptide[ww10], 5, 6), "XX", 
    substr(train.data$peptide[ww10], 7, 10))
  
  train.data$peprep[ww11] = paste0(
    substr(train.data$peptide[ww11], 1, 4), "XX", 
    substr(train.data$peptide[ww11], 5, 7), "XX", 
    substr(train.data$peptide[ww11], 8, 11))
  
  train.data$peprep[ww12] = paste0(
    substr(train.data$peptide[ww12], 1, 4), "XX", 
    substr(train.data$peptide[ww12], 5, 8), "X", 
    substr(train.data$peptide[ww12], 9, 12))
  
  train.data$peprep[ww13] = paste0(
    substr(train.data$peptide[ww13], 1, 4), "X", 
    substr(train.data$peptide[ww13], 5, 9), "X", 
    substr(train.data$peptide[ww13], 10, 13))
  
  train.data$peprep[ww14] = paste0(
    substr(train.data$peptide[ww14], 1, 4),  "X", 
    substr(train.data$peptide[ww14], 5, 14))
  
  train.data$peprep[ww15] = train.data$peptide[ww15]

  sum(is.na(train.data$peprep))
  
  # Create data that is 8-11 peptides and 12-15 peptides only 
  train.data8  = train.data[which(train.data$length>=8 & 
                                    train.data$length<=11),]
  train.data12 = train.data[which(train.data$length >= 12 & 
                                    train.data$length <=15)]
  
  print("ratio of nonbinders:binders in 8-11")
  print(round(table(train.data8$binder)[1]/table(train.data8$binder)[2],2))
  
  print(table(train.data8$binder))
  print("ratio of nonbinders:binders in 12-15")
  print(round(table(train.data12$binder)[1]/table(train.data12$binder)[2],2))
  print(table(train.data12$binder))
  
  # Output data files 
  dir1= "../../data/SA_data/by_split"
  
  cols2kp = c("peprep", "binder", "HLA")
  colsnms = c("peptide", "binder", "HLA")
  
  fnm = sprintf("%s/train_v4_el_single_HLA_%d.txt", dir1, i)
  write.table(train.data[,..cols2kp], fnm, sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = colsnms)
  system(sprintf("gzip %s", fnm))
  
  fnm = sprintf("%s/train_v4_el_8to11_single_HLA_%d.txt", dir1, i)
  write.table(train.data8[,..cols2kp], fnm, sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = colsnms)
  system(sprintf("gzip %s", fnm))
  
  fnm = sprintf("%s/train_v4_el_12to15_single_HLA_%d.txt", dir1, i)
  write.table(train.data12[,..cols2kp], fnm, sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = colsnms)
  system(sprintf("gzip %s", fnm))
  
  # Create dataset without any duplications (same peptide binding to 
  # at least 2 different HLA)
  # check for any duplicated peptides, remove all instances INCLUDING original
  print("Number of duplicated peptides")
  ndup = sum(duplicated(train.data$peptide) | 
               duplicated(train.data$peptide, fromLast=TRUE))
  print(ndup)
  temp = train.data[duplicated(train.data$peptide) | 
                      duplicated(train.data$peptide, fromLast=TRUE),]
  print(dim(temp))
  #b = as.data.frame(table(temp$peptide, temp$HLA))
  #c = b[which(b$Freq>0), ]
  #e = c[order(c$Var1),]
  #write.table(e, "TCR_HLA_table_train.txt", sep = " ", 
  #            quote= FALSE, row.names = TRUE, col.names = TRUE)
  #system("gzip TCR_HLA_table_train.txt")
  #table(table(e$Var1))

}

# --------------------------------------------------------------------------
# Create Final Training and Validation file
# --------------------------------------------------------------------------
train_fin = NULL
# create SA training data 
dir1= "../../data/SA_data/by_split"
dir2= "../../data/SA_data"

#train data
for(i in 0:3){
  temp = fread(sprintf("%s/train_v4_el_single_HLA_%d.txt.gz", dir1, i))
  train_fin = rbind(train_fin, temp)
}
fnm = sprintf("%s/train_v4_el_single_HLA.txt", dir2)
write.table(train_fin, fnm, sep="\t", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
system(sprintf("gzip %s", fnm))

#validation data
temp4 = fread(sprintf("%s/train_v4_el_single_HLA_4.txt.gz", dir1))
fnm = sprintf("%s/validate_v4_el_single_HLA.txt", dir2)

write.table(temp4, fnm, sep="\t", quote = FALSE, 
            row.names = FALSE, col.names = colsnms)
system(sprintf("gzip %s", fnm))

train_fin = fread(sprintf("%s/train_v4_el_single_HLA.txt.gz", dir2))
table(train_fin$binder)
length(unique(train_fin$HLA))

validate_fin = fread(sprintf("%s/validate_v4_el_single_HLA.txt.gz", dir2))
table(validate_fin$binder) 
length(unique(validate_fin$HLA))

table(train_fin$binder) + table(validate_fin$binder)

# --------------------------------------------------------------------------
# format MHC_pseudo.dat
# --------------------------------------------------------------------------
mhc = read.delim(paste0(dir0, "MHC_pseudo.dat"), sep = "", header = FALSE)
head(mhc)

write.table(mhc, "../../data/NetMHCpan4_1_train/MHC_pseudo.txt", sep = " ", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# --------------------------------------------------------------------------
# format allelist
# --------------------------------------------------------------------------
allist = read.delim(paste0(dir0, "allelelist"), sep = "", header = FALSE)
dim(allist)
head(allist)

write.table(allist, "../../data/NetMHCpan4_1_train/allelelist.txt", sep = " ", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

sessionInfo()
q(save="no")

