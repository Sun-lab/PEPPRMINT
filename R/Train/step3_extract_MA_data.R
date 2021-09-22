#---------------------------------------------------------------------------------
# Step 4. extract MA data for MixPep (NetMHCpan-4.1)
# Note: For MA HLA-I EL peptide data 
#---------------------------------------------------------------------------------
library(data.table)
library(stringi)
library(stringr)
library(ggplot2)
theme_set(theme_classic())

dir0 = "../../data/NetMHCpan4_1_train"

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

#table of cell_lines
hla.list = strsplit(allele.list$hla, split=",")
n.hla = sapply(hla.list, length)
table(n.hla)
table(n.hla > 1)

cell.line = allele.list$cell_line[n.hla > 1]

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
dat0[1:2,]

colnames(dat0) = c("peptide", "binder", "HLA")

# 1) remove SA lines
alleles0 = unique(dat0$HLA)
head(alleles0)

alleles1 = alleles0[grep("HLA-", alleles0, fixed=TRUE)]
length(alleles1)
alleles2 = alleles0[grep("H-2", alleles0, fixed=TRUE)]
alleles2
alleles3 = alleles0[grep("DLA", alleles0, fixed=TRUE)]
alleles3
alleles4 = alleles0[grep("Mamu", alleles0, fixed=TRUE)]
alleles4

#remove cell-lines
alleles5 = c("A10", "A11-A11", "A12-A15", "A14", "A15-A15", "A18", 
             "A19-A19", "A20-A20", "EBL")


alleles_sup3 = c(alleles1, alleles2, alleles3, alleles4, alleles5)
length(alleles_sup3)

cell.line = allele.list$cell_line[which(!allele.list$cell_line %in% alleles_sup3)]
length(cell.line)

dat1 = dat0[which(dat0$HLA %in% cell.line),]
length(unique(dat1$HLA))
table(dat1$binder)

allele.list1 = allele.list[which(allele.list$cell_line %in% dat1$HLA),]
dim(allele.list1)

hla.list1 = strsplit(allele.list1$hla, split=",")
n.hla = sapply(hla.list1, length)
table(n.hla)

table(dat0$HLA[which(dat0$HLA %in% cell.line)])
length(unique(dat0$HLA[which(dat0$HLA %in% cell.line)]))


df1 = data.frame(n_sample = as.numeric(table(dat1$HLA)))
g1 = ggplot(df1, aes(x=log10(n_sample))) + xlab("log10(# of peptides per sample)") + 
  geom_histogram(color="darkblue", fill="lightblue", bins=20)
pdf("figures/hist_n_sample_MA.pdf", width=3, height=2)
g1
dev.off()

summary(df1$n_sample)

# --------------------------------------------------------------------------
# create MA training dataset
# --------------------------------------------------------------------------
set.seed(1999)

cell.line = allele.list$cell_line[which(!allele.list$cell_line %in% alleles_sup3)]
length(cell.line)

#cell.lineMA = allele.list$cell_line[n.hla > 1]

for(i in 0:4){
  train.data0 = fread(sprintf("%s/c00%d_el", dir0, i))

  names(train.data0) = c("peptide", "binder", "cell_line")
  #print(dim(train.data0))
  if(nrow(train.data0)==0){
    print("ERROR; NO DATA in given file")
    break
  }
  
  #subset data to cell lines in allele.list
  train.data = train.data0[which(train.data0$cell_line %in% cell.line),]
  #print(dim(train.data))
  if(nrow(train.data)==0){
    print("ERROR: NO data in given cell lines")
    break
  }
  
  table(train.data$cell_line)
  
  table(train.data$binder)
  
  print(sprintf("i=%d", i))
  print(dim(train.data))
  
  #Extrat peptides of length 8-15
  ww1 = train.data$cell_line %in% cell.line & stri_length(train.data$peptide)<=15 &
    stri_length(train.data$peptide)>=8
  
  if(sum(ww1)==0){
    print("ERROR: NO Peptides of length 8 - 15")
    break
  }
  
  train.data = train.data[which(ww1),]
  train.data$length = nchar(train.data$peptide)
  
  print('MA data distribution before duplication filtering')
  print(table(train.data$binder))
  print(table(train.data$cell_line))
  print(table(train.data$length))
  
  #check for duplication for peptide+HLA, remove all instances except original in MA data
  print(dim(train.data))
  if(sum(duplicated(train.data[, c("peptide", "cell_line")]))>0){
    temp1 = train.data[, c("peptide", "cell_line")]
    train.data = train.data[!duplicated(temp1),]
    print(paste("Removed", sum(duplicated(temp1)),
                "duplicated peptide+cell_line paired observations in training. 
                Note: may need to confirm peptide+cell_line pair is not binder and nonbinder", sep = " "))
    # REMOVE PEPTIDES Considered binder and nonbinder!!!!!!!
  }else{
    print("Note: No duplicated peptide+cell_line paired observations in training" )
  }
  
  # DO NOT NEED TO check for any duplicated peptides binding to multiple cell lines
  
  #print(table(train.data$binder))
  #print(length(unique(train.data$cell_line)))
  #print(table(train.data$cell_line))
  #print(table(train.data$binder, train.data$cell_line))

  print('MA data distribution after duplication filtering')
  print(dim(train.data))
  print(table(train.data$binder))
  print(table(train.data$cell_line))
  print(table(train.data$pep_length))
  
  
  #create 15 length representation for training data 
  #print(table(train.data$cell_line))
  print(table(train.data$length))
  train.data$peprep = NA
  train.data$peprep[which(train.data$length==8)] = paste(
    substr(train.data$peptide[which(train.data$length==8)], 1, 4), "XXXXXXX", 
    substr(train.data$peptide[which(train.data$length==8)], 5, 8), sep="")
  
  train.data$peprep[which(train.data$length==9)] = paste(
    substr(train.data$peptide[which(train.data$length==9)], 1, 4), "XXX", 
    substr(train.data$peptide[which(train.data$length==9)], 5, 5), 
    "XXX", substr(train.data$peptide[which(train.data$length==9)], 6, 9),
    sep="")
  
  train.data$peprep[which(train.data$length==10)] = paste(
    substr(train.data$peptide[which(train.data$length==10)], 1, 4), "XXX",
    substr(train.data$peptide[which(train.data$length==10)], 5, 6), "XX",
    substr(train.data$peptide[which(train.data$length==10)], 7, 10), sep="")
  
  train.data$peprep[which(train.data$length==11)] = paste(
    substr(train.data$peptide[which(train.data$length==11)], 1, 4), "XX", 
    substr(train.data$peptide[which(train.data$length==11)], 5, 7), "XX",
    substr(train.data$peptide[which(train.data$length==11)], 8, 11), sep="")
  
  train.data$peprep[which(train.data$length==12)]=paste(
    substr(train.data$peptide[which(train.data$length==12)], 1, 4), "XX", 
    substr(train.data$peptide[which(train.data$length==12)], 5, 8), "X",
    substr(train.data$peptide[which(train.data$length==12)], 9, 12), sep="")
  
  train.data$peprep[which(train.data$length==13)]=paste(
    substr(train.data$peptide[which(train.data$length==13)], 1, 4), "X", 
    substr(train.data$peptide[which(train.data$length==13)], 5, 9), "X",
    substr(train.data$peptide[which(train.data$length==13)], 10, 13), sep="")
  
  train.data$peprep[which(train.data$length==14)] = paste(
    substr(train.data$peptide[which(nchar(train.data$peptide)==14)], 1, 4), "X", 
    substr(train.data$peptide[which(nchar(train.data$peptide)==14)], 5, 14), sep="")
  
  train.data$peprep[which(train.data$length==15)]=
    train.data$peptide[which(train.data$length==15)]
  
  sum(is.na(train.data$peprep))
  
  print(table(train.data$binder))
  print('Proportion of nonbinders to binders in MA data')
  print(table(train.data$binder)[1]/table(train.data$binder)[2])
  head(train.data)
  
  # Output data files 
  dir1= "../../data/MA_data"
  
  cols2kp = c("peprep", "binder", "cell_line")
  colsnms = c("peptide", "binder", "cell_line")
  
  fnm = sprintf("%s/train_v4_el_multi_only_HLA_%d.txt", dir1, i)
  write.table(train.data[,..cols2kp], fnm, sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = colsnms)
  system(sprintf("gzip %s", fnm))
  
  # add SA training data
  SA = fread(sprintf("../../data/SA_data/by_split/train_v4_el_single_HLA_%d.txt.gz", i))
  
  head(SA)
  colnames(SA) = c("peprep", "binder", "cell_line")
  SA$length = nchar(str_replace_all(SA$peprep, "X", ""))
  print("SA training data")
  print(dim(SA))
  print(table(SA$binder))
  print(length(unique(SA$cell_line)))
  print(table(SA$length))
  
  train.datasub = train.data[,c("peprep", "binder", "cell_line", "length")]
  
  train.data.fin = rbind(train.datasub, SA)
  head(train.data.fin)
  
  print("Final data characteristics (MA + SA) ")
  print(dim(train.data.fin))
  print(table(train.data.fin$binder))
  print("Number of cell lines")
  print(length(unique(train.data.fin$cell_line)))
  print('Proportion of nonbinders to binders in MA data')
  print(table(train.data.fin$binder)[1]/table(train.data.fin$binder)[2])
  print(table(train.data.fin$length))
  
  #subset to 8-11 and 12-15
  ds_8 = train.data.fin[which(train.data.fin$length>=8 & train.data.fin$length<=11),]
  print("8-11 peptides only: # cell lines,")
  dim(ds_8)
  print(length(unique(ds_8$cell_line)))
  print(table(ds_8$length))
  
  ds_12 = train.data.fin[which(train.data.fin$length>=12 ),]
  print("12-15 peptides only: # cell lines,")
  dim(ds_12)
  print(length(unique(ds_12$cell_line)))
  print(table(ds_12$length))
  
  
  
  # Output data files 
  dir1= "../../data/MA_data"
  
  cols2kp = c("peprep", "binder", "cell_line")
  colsnms = c("peptide", "binder", "cell_line")
  
  fnm = sprintf("%s/train_v4_el_multi_HLA_%d.txt", dir1, i)
  write.table(train.data.fin[,..cols2kp], fnm, sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = colsnms)
  system(sprintf("gzip %s", fnm))
  
  fnm = sprintf("%s/train_v4_el_8to11_multi_HLA_%d.txt", dir1, i)
  write.table(ds_8[,..cols2kp], fnm, sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = colsnms)
  system(sprintf("gzip %s", fnm))
  
  fnm = sprintf("%s/train_v4_el_12to15_multi_HLA_%d.txt", dir1, i)
  write.table(ds_8[,..cols2kp], fnm, sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = colsnms)
  system(sprintf("gzip %s", fnm))
  
}

sessionInfo()
q(save="no")

