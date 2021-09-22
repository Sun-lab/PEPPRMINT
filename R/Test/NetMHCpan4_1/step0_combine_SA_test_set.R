# ------------------------------------------------------------------------------------
# test dataset for HLA-I 
# ------------------------------------------------------------------------------------

dir0 = "../../../data/test_data/NetMHCpan4_1/"
library(stringr)
library(data.table)

fall = list.files("../../../data/test_data/NetMHCpan4_1/",
                  pattern=".txt")
fall

test = NULL
for(f1 in fall){
  print(f1)
  temp = fread(paste(dir0, f1, sep =""))
  print(dim(temp))
  
  test = rbind(test, temp)
}

dim(test)
colnames(test) = c("peptide", "binder", "HLA")
head(test)

table(test$binder)
table(test$HLA)

#Put in 15-length representation
#create 15 length representation for training data 
makerep = function(train.data){
  #create 15 length representation for training data 
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

test_pad = makerep(test) 
head(test_pad)

# ------------------------------------------------------------------------------------
# Output data files 
# ------------------------------------------------------------------------------------

cols2kp = c("peprep", "binder", "HLA")
colsnms = c("peptide", "binder", "HLA")

fnm = sprintf("%s/NetMHCpan41_test_el_single_HLA.txt", dir0)
write.table(test_pad[,..cols2kp], fnm, sep="\t", quote = FALSE, 
            row.names = FALSE, col.names = colsnms)
system(sprintf("gzip %s", fnm))


sessionInfo()
q(save="no")

