# ---------------------------------------------------------------------------
# Figures 
# ---------------------------------------------------------------------------

#Import hla specific resutls 
library(data.table)
hla_auc_results = fread("../../../results/test_data/NetMHCpan4_1/by_hla_auc_PEPPRMINT.txt")
hla_auc_results_net = fread("../../../results/test_data/NetMHCpan4_1/by_hla_auc_NetMHCpan4_1.txt")

hla_auc_results_comb = merge(hla_auc_results, hla_auc_results_net, by = "hla")
dim(hla_auc_results_comb)

#add MHCflurry-2.0 results (Supplementary Table 8 in NetMHCpan-4.1)
hla_auc_results_comb$AUC_fl = c(0.97215,
                                0.93177,
                                0.97386,
                                0.96266,
                                0.93949,
                                0.87331,
                                0.92271,
                                0.92205,
                                0.91348,
                                0.91974,
                                0.94563,
                                0.96796,
                                0.92158,
                                0.89455,
                                0.87285,
                                0.93989,
                                0.93493,
                                0.91402,
                                0.93979,
                                0.96589,
                                0.97423,
                                0.95250,
                                0.92374,
                                0.88834,
                                0.92332,
                                0.97796,
                                0.95931,
                                0.96262,
                                0.92001,
                                0.97703,
                                0.94191,
                                0.80148,
                                0.90822,
                                0.95192,
                                0.94975,
                                0.95953)
# ---------------------------------------------------------------------------
# Scatter plot by HLA
# ---------------------------------------------------------------------------
library(ggrepel)

# Comparison to NetMHCpan4.1

bp = ggplot(hla_auc_results_comb, aes(x = AUC_net, 
                                     y = AUC, color = hla)) + 
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
        #legend.key.size = unit(0.5,"line"),
        legend.key.height= unit(0.5, 'cm'), 
        legend.key.width= unit(2.5, 'cm'),
        axis.ticks.length = unit(0.6, "cm"))+
  guides(color = guide_legend(nrow = 6))     


# Comparison to MHCflurry-2.0
cp = ggplot(hla_auc_results_comb, aes(x = AUC_fl, 
                                     y = AUC, color = hla)) + 
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
        legend.text = element_text(size = 40), 
        legend.title =element_blank(),
        #legend.key.size = unit(2,"line"),
        legend.key.height= unit(0.5, 'cm'), 
        legend.key.width= unit(4, 'cm'),
        axis.ticks.length = unit(0.6, "cm"))+
  guides(color = guide_legend(nrow = 4))     

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
comb_p = grid.arrange(bp, cp, ncol=2, nrow = 1, 
                      layout_matrix = rbind(c(1,2), c(3,3)),
                      widths = c(3, 3), heights = c(3))

l = grid.arrange(legend, ncol=1, nrow = 1,
                 widths = c(20), heights = c(4))

# 5. Save as png
ggsave(
  "../../../figures/PEPPRMINT/hla_AUC_compare.pdf",
  plot = comb_p,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 46,
  height = 22,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

ggsave(
  "../../../figures/PEPPRMINT/hla_legend.pdf",
  plot = l,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 26,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

#ROC Net vs. PEPPRMINT
pdf(file="../../figures/PEPPRMINT/AUC_compare_Sep3.pdf",
    width=5, height=5, units = 'in', res = 300)
par(mar=c(8, 5, 5, 5), xpd=TRUE)
plot(netroc, color = "darkblue", auc.main=FALSE,  lty = 3, bty= "n", 
     cex.axis=1, cex.lab = 1, main = "", cex.main = 2)
plot(roc1, color = "darkorange", add=TRUE)
legend(cex = .8,"bottom", inset=c(0,-.5), c(paste0("PEPPRMINT (AUC = ", round(roc1$auc,3), 
                                                   ", PPV = ", round(ppv,3), ")"),
                                            paste0("NetMHCpan-4.1 (AUC = ", round(netroc$auc,3), 
                                                   ", PPV = ", round(ppvnet,3), ")")
), 
col=c("darkorange", "darkblue"), lwd=c(2,2), lty = c(1,2), bty="n")

# SUPPLEMENTARY TABLE 
colnames(hla_auc_results_comb) = c("HLA", "AUC_PEPPRMINT", "PPV_PEPPRMINT", 
                                   "AUC_NetMHCpan4.1", "PPV_NetMHCpan4.1", 
                                   "AUC_MHCflurry2.0")
hla_auc_results_comb$AUC_PEPPRMINT= round(hla_auc_results_comb$AUC_PEPPRMINT, 3)
hla_auc_results_comb$AUC_NetMHCpan4.1= round(hla_auc_results_comb$AUC_NetMHCpan4.1, 3)
hla_auc_results_comb$AUC_MHCflurry2.0= round(hla_auc_results_comb$AUC_MHCflurry2.0, 3)

hla_auc_results_comb$PPV_PEPPRMINT= round(hla_auc_results_comb$PPV_PEPPRMINT, 3)
hla_auc_results_comb$PPV_NetMHCpan4.1= round(hla_auc_results_comb$PPV_NetMHCpan4.1, 3)

pdf("../../../results/test_data/SuppTable2_SAtest_performance_by_sample_method_comparison.pdf",
    height=11, width=10)
grid.table(hla_auc_results_comb, rows= NULL)
dev.off()

dev.off()


