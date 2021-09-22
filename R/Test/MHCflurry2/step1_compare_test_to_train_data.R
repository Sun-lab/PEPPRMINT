# ------------------------------------------------------------------------------------
# Comparison with training data used for PEPPRMINT 
#   Note: PEPPRMINT training data is the human allele subset of the NetMHCpan-4.1 data
# ------------------------------------------------------------------------------------
dir0 = "../../../data/MA_data/"
fall = list.files(dir0, pattern="train_v4_el_multi_HLA_")
fall

# ------------------------------------------------------------------------------------
# Compare with NetMHCpan-4.1 training data 
# ------------------------------------------------------------------------------------
nettrain = fread(paste0(dir0, fall[[1]]))
fall1 = fall[-1]

for(f in fall1){
  print(f)
  split = fread(paste0(dir0, f))
  nettrain = rbind(nettrain,split)
}

dim(nettrain)
head(nettrain)
nettrain$Peptide = str_replace_all(nettrain$peptide, "X", "")
head(nettrain)

net1 = nettrain[which(nettrain$binder==1),]
dim(net1)

net1sub = net1[,c("peptide", "Peptide")]
net1sub = distinct(net1sub)
dim(net1sub)

# MHCflurry-2.0 MULTIALLELIC-RECENT dataset 
mhcflurry= fread("../../../data/test_data/MHCflurry2/MHCflurry2_test_el_multi_HLA.txt.gz")
flurry_bind = mhcflurry[mhcflurry$binder==1,]
dim(flurry_bind)

over = merge(flurry_bind, net1sub, by = c("peptide", "Peptide"))
head(over)
dim(over)

# output set of peptides that are not in NetMHCpan-4.1 training data
over_pep  = over$Peptide
a = flurry_bind$Peptide %in% over_pep 
dis_pep = flurry_bind[!a,]

dis_pep$sample = paste(dis_pep$cell_line,
                                 dis_pep$peptide, sep = ";")
head(dis_pep)

dis_pep$pep_line = paste(dis_pep$cell_line,
                       dis_pep$Peptide, sep = ";")
head(dis_pep)



#save 
fnm = "../../../data/test_data/MHCflurry2/MHCflurry2_test_el_multi_HLA_no_overlap.txt"
write.table(dis_pep[,c("sample", "binder", "pep_line")], fnm,
            sep = "\t", quote = FALSE, row.names = FALSE, 
            col.names = c("peptide", "y_true", "pep_line"))



over1 = merge(flurry_bind, net1, by = c("peptide", "Peptide"))
head(over1)
dim(over1)
table(over1$cell_line.x)
table(over1$cell_line.x, over1$cell_line.y)

# conclusion: half of peptides are in test set (15150 of 27007), however they do not match 1-1 in cell-line

sessionInfo()
q(save="no")

