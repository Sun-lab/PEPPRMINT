# ---------------------------------------------------------------------------
# Figures for MA test data 
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Import PEPPRMINT Results 
# ---------------------------------------------------------------------------
all_max1 = fread("../../../results/test_data/MHCflurry2/PEPPRMINT_MHCflurry2_pred_Aug5.txt")

#ROC
roc1 = roc.curve(scores.class0 = all_max1$cv_score, 
                 weights.class0 = all_max1$y_true, 
                 curve=TRUE)
roc1
#plot(roc1)

#PPV
temp1 <- all_max1[with(all_max1,order(-cv_score)),]
hits = table(all_max1$y_true)[2][[1]]
temp2 = temp1[1:hits,]
dim(temp2)

ppv = sum(temp2$y_true)/hits
ppv

cl_auc_results = fread("../../../results/test_data/MHCflurry2/cell_line_performance_PEPPRMINT.txt")
dim(cl_auc_results)
head(cl_auc_results)

# ---------------------------------------------------------------------------
# Import NetMHCpan4-1 Results 
# ---------------------------------------------------------------------------
net_pred = fread("../../../results/test_data/MHCflurry2/NetMHCpan4_1_MHCflurry2_test_pred.txt")

#AUC
netroc = roc.curve(scores.class0 = net_pred$max_pred, 
                   weights.class0 = net_pred$binder, curve = TRUE)

netroc

# PPV
hits = table(net_pred$binder)[2][[1]]
temp1 <- net_pred[with(net_pred,order(-max_pred)),]
temp2 = temp1[1:hits,]
dim(temp2)

ppvnet = sum(temp2$binder)/hits
ppvnet

cl_net_auc =  read.delim("../../../results/test_data/MHCflurry2/NetMHCpan4_1_MHCflurry2_cell_line_performance.txt", 
                                sep = "\t")

cl_net_auc$sample_id = as.character(cl_net_auc$sample_id)
cl_net_auc$NetMHCpan4.1_AUC = as.numeric(as.character(cl_net_auc$NetMHCpan4.1_AUC))
cl_net_auc$NetMHCpan4.1_PPV = as.numeric(as.character(cl_net_auc$NetMHCpan4.1_PPV))

cl_net_auc$sample_id[is.na(cl_net_auc$sample_id)] = "29_14-TISSUE"
cl_net_auc$sample_id[which(cl_net_auc$sample_id=="29_14-TISSUE")] = "29/14-TISSUE"
cl_net_auc$sample_id[which(cl_net_auc$sample_id=="637_13-TISSUE")] = "637/13-TISSUE"

cl_auc_results_comb = merge(cl_net_auc, 
                            cl_auc_results, by = intersect(colnames(cl_auc_results), 
                                                           colnames(cl_net_auc)))
dim(cl_auc_results_comb)
head(cl_auc_results_comb)

table(cl_auc_results_comb$AUC>cl_auc_results_comb$NetMHCpan4.1_AUC)

# ---------------------------------------------------------------------------
# Import MHCflurry-2.0
# ---------------------------------------------------------------------------
fl_res= read.delim("../../../data/test_data/MHCflurry2/performance_results.txt", 
                   sep = "\t")


head(fl_res )

#add MHCflurry2.0 results 
cl_auc_results_comb = merge(fl_res[,c("sample_id", "AUC..MHCflurry.2.0.AP_plus_flanks",
                                      "PPV..MHCflurry.2.0.AP_plus_flanks")], 
                            cl_auc_results_comb, 
                            by = intersect(colnames(cl_auc_results_comb), 
                                           colnames(fl_res)))

dim(cl_auc_results_comb)
head(cl_auc_results_comb)


# ---------------------------------------------------------------------------
# Figures 
# ---------------------------------------------------------------------------
library(ggrepel)
# Comparison to NetMHCpan4.1

bp = ggplot(cl_auc_results_comb, aes(x = NetMHCpan4.1_AUC, 
                                     y = AUC, color = sample_id)) + 
  geom_point(size = 20) + 
  #ggtitle("AUC Comparison to NetMHCpan-4.1" ) + 
  scale_x_continuous(name="NetMHCpan-4.1", limits=c(.75, 1)) +
  scale_y_continuous(name="PEPPRMINT", limits=c(.75, 1)) + 
  geom_segment(aes(x = 0.75, y = 0.75, xend = 1, yend = 1),
               colour='black', lwd = 2) + 
  theme_classic() + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = 70), 
        axis.title = element_text(size = 85), 
        plot.title = element_text(size = 30, face ="bold"), 
        legend.text = element_text(size = 40), 
        legend.title = element_blank(), 
        #legend.key.size = unit(3,"cm"),
        legend.key.height= unit(0.5, 'cm'), 
        legend.key.width= unit(4, 'cm'),
        axis.ticks.length = unit(0.6, "cm"))+
  guides(color = guide_legend(nrow = 5))     


# Comparison to MHCflurry-2.0
cp = ggplot(cl_auc_results_comb, aes(x = AUC..MHCflurry.2.0.AP_plus_flanks, 
                                     y = AUC, color = sample_id)) + 
  geom_point(size = 20) + 
  #ggtitle("AUC Comparison to MHCflurry-2.0" ) + 
  scale_x_continuous(name="MHCflurry-2.0", limits=c(.75, 1)) +
  scale_y_continuous(name="PEPPRMINT", limits=c(.75, 1)) + 
  geom_segment(aes(x = 0.75, y = 0.75, xend = 1, yend = 1),
               colour='black', lwd = 2) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 70), 
        axis.title = element_text(size = 85), 
        plot.title = element_text(size = 30, face ="bold"), 
        legend.text = element_text(size = 22), 
        legend.title = element_blank(), 
        legend.key.height= unit(0.5, 'cm'), 
        legend.key.width= unit(4, 'cm'),
        axis.ticks.length = unit(0.6, "cm"))+
  guides(color = guide_legend(nrow = 5))     

# 2. Save the legend
#+++++++++++++++++++++++
library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(bp)

# 3. Remove the legend from the b plot
#+++++++++++++++++++++++
bp <- bp + theme(legend.position="none")
cp <- cp + theme(legend.position="none")
# 4. Arrange ggplot2 graphs with a specific width
#comb_p = grid.arrange(bp, cp, legend, ncol=2, nrow = 2, 
#                      layout_matrix = rbind(c(1,2), c(3,3)),
#                      widths = c(5, 5), heights = c(4, 1))

comb_p1 = grid.arrange(bp, cp, ncol=2, nrow = 1, 
                      layout_matrix = rbind(c(1,2), c(3,3)),
                      widths = c(5, 5), heights = c(4))
l = grid.arrange( legend, ncol=1, nrow = 1,
                 widths = c(20), heights = c(4))

# 5. Save as png
ggsave(
  "../../../figures/PEPPRMINT/cell_line_AUC_compare_nolegend.pdf",
  plot = comb_p1,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 44,
  height = 22,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

ggsave(
  "../../../figures/PEPPRMINT/legend.pdf",
  plot = l,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 28,
  height = 4,
  dpi = 300,
  limitsize = FALSE,
  bg = NULL
)



#ROC Net vs. PEPPRMINT
pdf(file="../../figures/PEPPRMINT/AUC_compare_Sep3.pdf",
    width=5, height=5)
par(mar=c(8, 5, 5, 5), xpd=TRUE)
plot(netroc, color = "darkblue", auc.main=FALSE,  lty = 2, bty= "n", 
     cex.axis=1, cex.lab = 1, main = "", cex.main = 2)
plot(roc1, color = "darkorange", add=TRUE)
legend(cex = .8,"bottom", inset=c(0,-.5), c(paste0("PEPPRMINT (AUC = ", round(roc1$auc,3), 
                                                   ", PPV = ", round(ppv,3), ")"),
                                            paste0("NetMHCpan-4.1 (AUC = ", round(netroc$auc,3), 
                                                   ", PPV = ", round(ppvnet,3), ")")
), 
col=c("darkorange", "darkblue"), lwd=c(2,2), lty = c(1,2), bty="n")

dev.off()

# ---------------------------------------------------------------------------
# ---- SUPPLEMENTARY TABLES 
# ---------------------------------------------------------------------------
colnames(cl_auc_results_comb) = c("sample", "AUC_MHCflurry2", "PPV_MHCflurry2",
                                  "AUC_NetMHCpan4.1", "PPV_NetMHCpan4.1", 
                                  "AUC_PEPPRMINT", "PPV_PEPPRMINT")
cl_auc_results_comb$AUC_PEPPRMINT= round(cl_auc_results_comb$AUC_PEPPRMINT, 3)
cl_auc_results_comb$AUC_NetMHCpan4.1= round(cl_auc_results_comb$AUC_NetMHCpan4.1, 3)
cl_auc_results_comb$AUC_MHCflurry2= round(cl_auc_results_comb$AUC_MHCflurry2, 3)

cl_auc_results_comb$PPV_PEPPRMINT= round(cl_auc_results_comb$PPV_PEPPRMINT, 3)
cl_auc_results_comb$PPV_NetMHCpan4.1= round(cl_auc_results_comb$PPV_NetMHCpan4.1, 3)
cl_auc_results_comb$PPV_MHCflurry2= round(cl_auc_results_comb$PPV_MHCflurry2, 3)

cl_auc_results_comb = cl_auc_results_comb[,c("sample", "AUC_PEPPRMINT", "AUC_NetMHCpan4.1",
                                             "AUC_MHCflurry2", "PPV_PEPPRMINT",
                                             "PPV_NetMHCpan4.1", "PPV_MHCflurry2")]

pdf("../../../results/test_data/SuppTable1_MAtest_performance_by_sample_method_comparison.pdf",
    height=6, width=13)
grid.table(cl_auc_results_comb, rows= NULL)
dev.off()

write.table(cl_auc_results_comb,
            , 
            sep = " ", quote = FALSE, col.names = TRUE, row.names = FALSE)
sessionInfo()
q(save="no")



########################################## ---------------------------------------------------------------------------
# Supplementary PPV plot 
ggplot(cl_auc_results_comb, aes(x = NetMHCpan4.1_PPV, 
                                y = PPV, color = sample_id)) + 
  geom_point(size = 3) + 
  ggtitle("PPV Cell-line Comparison to NetMHCpan-4.1") + 
  scale_x_continuous(name="NetMHCpan-4.1", limits=c(0, .7)) +
  scale_y_continuous(name="PEPPRMINT", limits=c(0,.7)) + 
  geom_segment(aes(x = 0, y = 0, xend = .7, yend = .7),
               colour='black',  lwd = 0.5) #+ 
#geom_text_repel()

ggplot(cl_auc_results_comb, aes(x = PPV..MHCflurry.2.0.AP_plus_flanks, 
                                y = PPV, color = sample_id)) + 
  geom_point(size = 3) + 
  ggtitle("PPV Cell-line Comparison to MHCflurry-2.0") + 
  scale_x_continuous(name="MHCflurry-2.0", limits=c(0, .7)) +
  scale_y_continuous(name="PEPPRMINT", limits=c(0,.7)) + 
  geom_segment(aes(x = 0, y = 0, xend = .7, yend = .7),
               colour='black', lwd =0.5) #+ 
#geom_text_repel()
