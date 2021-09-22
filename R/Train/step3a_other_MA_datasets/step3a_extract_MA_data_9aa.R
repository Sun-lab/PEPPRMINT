#---------------------------------------------------------------------------------
# Step 4. extract 9 AA MA data
# Note: For MA HLA-I EL peptide data 
#---------------------------------------------------------------------------------

library(data.table)
library(stringi)
library(stringr)

dir0 = "../../../data/NetMHCpan4_1_train"

# --------------------------------------------------------------------------
# Extract EL data that are associated with a multiple HLA alleles
# --------------------------------------------------------------------------
# this file has been modified to replace tab with space in 'step1'
allele.list = fread(sprintf("%s/allelelist.txt", dir0), header=FALSE)
dim(allele.list)
names(allele.list) = c("cell_line", "hla")
head(allele.list)

#remove duplciations 
allele.list[duplicated(allele.list$cell_line),]
allele.list[which(allele.list$cell_line=="BoLA-6:01301"),]
allele.list = allele.list[!duplicated(allele.list$cell_line),]
dim(allele.list)

# table of cell_lines
hla.list = strsplit(allele.list$hla, split=",")
n.hla    = sapply(hla.list, length)
table(n.hla)
table(n.hla > 1)

# cell lines with HLA alleles
hla_cell_line = grepl("HLA", allele.list$cell_line)
hla_in_hla    = grepl("HLA", allele.list$hla)
table(hla_cell_line, hla_in_hla, useNA="ifany")

# --------------------------------------------------------------------------
# Check MA data matches Supplementary Table 1 and 5
# --------------------------------------------------------------------------
# import all data
c0 = fread(paste0(dir0, "/c000_el"))
c1 = fread(paste0(dir0, "/c001_el"))
c2 = fread(paste0(dir0, "/c002_el"))
c3 = fread(paste0(dir0, "/c003_el"))
c4 = fread(paste0(dir0, "/c004_el"))

# combine all data 
dat0 = rbind(c0, c1)
dat0 = rbind(dat0, c2)
dat0 = rbind(dat0, c3)
dat0 = rbind(dat0, c4)
dim(dat0)
colnames(dat0) = c("peptide", "binder", "HLA")
dat0[1:2,]

cell_line_2use = allele.list$cell_line[hla_in_hla & (! hla_cell_line)]
length(cell_line_2use)

dat1 = dat0[which(dat0$HLA %in% cell_line_2use),]
dim(dat1)
dat1[1:2,]

length(unique(dat1$HLA))
table(dat1$binder)
length(unique(dat1$HLA))
sort(table(dat1$HLA))

# --------------------------------------------------------------------------
# create MA training dataset
# --------------------------------------------------------------------------
set.seed(1999)

for(i in 0:4){
  train.data0 = fread(sprintf("%s/c00%d_el", dir0, i))

  names(train.data0) = c("peptide", "binder", "cell_line")

  if(nrow(train.data0)==0){
    print("ERROR; NO DATA in given file")
    break
  }
  
  train.data = train.data0[cell_line %in% cell_line_2use,]

  if(nrow(train.data)==0){
    print("ERROR: NO data in given cell lines")
    break
  }
  
  print(sprintf("i=%d", i))
  print(dim(train.data))
  
  #Extrat peptides of length 8-15
  pep_len = nchar(train.data$peptide)
  ww1 = pep_len <=15 & pep_len >= 8
  
  if(sum(ww1)==0){
    print("ERROR: NO Peptides of length 8 - 15")
    break
  }
  
  train.data = train.data[which(ww1),]
  train.data$length = nchar(train.data$peptide)
  
  print('MA data distribution before duplication filtering')
  print(table(train.data$binder))
  print(table(train.data$length))
  
  # check for duplication for peptide+HLA
  print(dim(train.data))
  if(sum(duplicated(train.data[, c("peptide", "cell_line")]))>0){
    temp1 = train.data[, c("peptide", "cell_line")]
    train.data = train.data[!duplicated(temp1),]
    print(paste("Removed", sum(duplicated(temp1)),
                "duplicated peptide+cell_line paired observations"))
  }else{
    print("Note: No duplicated peptide+cell_line paired observations" )
  }
  
  print('MA data distribution after duplication filtering')
  print(dim(train.data))
  print(table(train.data$binder))
  print(table(train.data$length))
  
  ww8  = which(train.data$length==8)
  ww9  = which(train.data$length==9)
  
  data_all = NULL
  data_new = NULL
  
  data8 = train.data[ww8,]
  data8$pep_core  = gsub("^(.{3})(.*)$", "\\1X\\2", data8$peptide)
  data8$start_pos = 4
  data_new = rbind(data_new, data8)
  
  data8 = train.data[ww8,]
  data8$pep_core = gsub("^(.{4})(.*)$", "\\1X\\2", data8$peptide)
  data8$start_pos = 5
  data_new = rbind(data_new, data8)
  
  data8 = train.data[ww8,]
  data8$pep_core = gsub("^(.{5})(.*)$", "\\1X\\2", data8$peptide)
  data8$start_pos = 6
  data_new = rbind(data_new, data8)
  
  # for each non-binder, randomly choose one start position
  data_n1 = data_new[binder==1,]
  data_n0 = data_new[binder==0,]
  
  data_n0$key = paste(data_n0$peptide, data_n0$cell_line, sep=":")
  stopifnot(all(table(data_n0$key) == 3))
  
  dim(data_n0)
  data_n0
  data_n0 = data_n0[, .SD[sample(x = .N, size = 1)], by = key]
  dim(data_n0)
  data_n0
  stopifnot(all(table(data_n0$key) == 1))
  data_n0[, key:=NULL]
  
  data_n1$weight = 1/3
  data_n0$weight = 1
  
  data_all = rbind(data_all, data_n1)
  data_all = rbind(data_all, data_n0)
  
  data9 = train.data[ww9,]
  data9$pep_core = data9$peptide
  data9$start_pos = 1
  data9$weight = 1
  data_all = rbind(data_all, data9)
  
  for(len1 in 10:15){
    data_l = train.data[which(train.data$length==len1),]
    len2rm = len1 - 9
    
    data_new = NULL
    for(k in 4:(len1 - 3 - (len2rm-1))){
      data_l$pep_core = data_l$peptide
      str_sub(data_l$pep_core, start=k, end=k+len2rm-1) = ""
      data_l$start_pos = k
      data_new = rbind(data_new, data_l)
    }
    
    data_n1 = data_new[binder==1,]
    data_n0 = data_new[binder==0,]
    
    data_n0$key = paste(data_n0$peptide, data_n0$cell_line, sep=":")
    stopifnot(all(table(data_n0$key) == 4))
    
    dim(data_n0)
    data_n0
    data_n0 = data_n0[, .SD[sample(x = .N, size = 1)], by = key]
    dim(data_n0)
    data_n0
    stopifnot(all(table(data_n0$key) == 1))
    data_n0[, key:=NULL]
    
    data_n1$weight = 1/4
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

  # add SA training data
  dir2 = "../../../data/SA_data/by_split"
  fnm  = sprintf("%s/train_v4_el_single_HLA_9AA_%d.txt.gz", dir2, i)
  SA   = fread(fnm)
  
  head(SA)
  colnames(SA)[3] = ("cell_line")
  
  data_all = rbind(data_all, SA)
  head(data_all)
  
  print("Final data characteristics (MA + SA) ")
  print(dim(data_all))
  print(table(data_all$binder))
  print("Number of cell lines")
  print(length(unique(data_all$cell_line)))
  print('Proportion of nonbinders to binders in MA data')
  print(table(data_all$binder)[1]/table(data_all$binder)[2])
  print(table(data_all$length))
  
  # Output data files 
  dir1= "../../../data/mixPep_data"
  
  fnm = sprintf("%s/train_v4_el_multi_HLA_9AA_%d.txt", dir1, i)
  fwrite(data_all, fnm, sep="\t")
  system(sprintf("gzip %s", fnm))
}

gc()
sessionInfo()
q(save="no")

