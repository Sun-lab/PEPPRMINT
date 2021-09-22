# ------------------------------------------------------------------------------------
# MHCflurry-2.0 Test data for PEPPRMINT
# O'Donnell et al (2020)
#   1. Create test data for PEPPRMINT Prediction 
#   2. Format allele list file 
# ------------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(stringr)

# ------------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------------

flurry = fread("../../../data/test_data/MHCflurry2/Data_S1.csv")

head(flurry)
table(flurry$sample_group)
table(flurry$hit)

# ------------------------------------------------------------------------------------
# MA test set
# ------------------------------------------------------------------------------------
flurry1 = flurry[which(flurry$sample_group=="MULTIALLELIC-RECENT"),]
dim(flurry1)

table(flurry1$hit)
table(flurry1$hit, flurry1$sample_id)

flurry_MA = flurry1[,c("peptide", "hit", "sample_id")]
flurry_MA$length = nchar(flurry_MA$peptide)
table(flurry_MA$hit, flurry_MA$length)
dim(flurry_MA)
head(flurry_MA)

makerep = function(train.data){
  #create 15 length representation for data 
  train.data$length = nchar(train.data$peptide)
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
  
  print(dim(train.data))
  
  
  return(train.data)
}

test_pad = makerep(flurry_MA) 
head(test_pad)
table(test_pad$hit)

# ------------------------------------------------------------------------------------
# Output test data files 
# ------------------------------------------------------------------------------------
dir1= "../../../data/test_data/MHCflurry2/"

cols2kp = c("peprep", "hit", "sample_id")
colsnms = c("peptide", "binder", "cell_line")

fnm = sprintf("%s/MHCflurry2_test_el_multi_HLA.txt", dir1)
write.table(test_pad[,..cols2kp], fnm, sep="\t", quote = FALSE, 
            row.names = FALSE, col.names = colsnms)
system(sprintf("gzip %s", fnm))


# ------------------------------------------------------------------------------------
# Create allele-list file
# ------------------------------------------------------------------------------------
flurry2 = flurry1[,c("sample_id", "hla")]
head(flurry2)
flurry3 = distinct(flurry2)
flurry3

#replace all spaces with ,
flurry3$hlas = str_replace_all(flurry3$hla, " ", ",")
#replace all * with nothing
flurry3$hlas = str_replace_all(flurry3$hlas, "\\*", "")
flurry3

#output allele list
fnm = sprintf("%s/allelelist_MHCflurry2.txt", dir1)
write.table(flurry3[,c("sample_id", "hlas")], fnm, sep=" ", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)



sessionInfo()
q(save="no")









