rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rstatix)

### IMPORT
ct_psnv <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k7_pSNV_codon_table.csv")

### BOXPLOT k7 DISTRIBUTIONS by category 

dt_box <- ct_psnv
dt_box$mutation2 <- dt_box$mutation
dt_box$mutation2[dt_box$mutation != "CG>TG"] <- paste0("  ", dt_box$mutation[dt_box$mutation != "CG>TG"], "  ")

p_box1 <- ggplot(dt_box[mutation != "CG>TG"], mapping = aes(x=mutation2, y=k7_mu_rate)) +
  # geom_violin(width=1.1, fill="grey") +
  # geom_boxplot(width=0.8, outlier.shape = NA, fill="white") +
  geom_boxplot(width=0.8, outlier.shape = NA, aes(fill=mutation2)) +
  ylab("7-mer substitution rate") +
  xlab("") +
  scale_fill_manual(values = rev(c("#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442"))) +
  # ggtitle("") +
  ylim(0, 0.1) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        plot.margin = unit(c(0.1,0.2,0.2,0.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p_box1

p_box2 <- ggplot(dt_box[mutation == "CG>TG"], mapping = aes(x=mutation2, y=k7_mu_rate)) +
  # geom_violin(width=1, fill="grey") +
  # geom_boxplot(width=0.7, outlier.shape = NA, fill="white") +
  geom_boxplot(width=0.7, outlier.shape = NA, aes(fill=mutation2)) +
  ylab("") +
  xlab("") +
  ggtitle("(B)") +
  scale_fill_manual(values = "#E69F00") +
  # scale_fill_brewer(palette="Set1") +
  # ylim(0, 0.1) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        plot.margin = unit(c(0.2,0.2,0.1,0.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p_box2
p_box <- grid.arrange(p_box2, p_box1, nrow = 2, ncol=1, heights=c(2, 7))

### BOXPLOT k7 DISTRIBUTIONS by category 

tmp <- ct_psnv[Mutation_type %in% c("nonsense", "synonymous", "missense") & mu_group %in% c("A>", "C>")]
tmp$mu_group <- "(A|C)>"
dt_plot <- rbind(tmp, ct_psnv[Mutation_type %in% c("nonsense", "synonymous", "missense") & mu_group %in% c("A>", "C>")])

stat.test <- dt_plot %>%
  group_by(mu_group) %>%
  wilcox_test(k7_mu_rate ~ Mutation_type, ref.group = "synonymous") 

stat.test <- stat.test %>%
  add_xy_position(x = "mu_group", dodge = 0.8)
stat.test <- as.data.table(stat.test)
stat.test$y.position <- c(0.09, (0.09+0.02),
                          0.06, (0.06+0.02),
                          0.126, (0.126+0.02))

p_wilcox <- ggplot(dt_plot, mapping = aes(x=mu_group, y=k7_mu_rate)) +
  geom_boxplot(width=0.6, outlier.shape = NA, position=position_dodge(0.8), aes(fill=Mutation_type)) +
  scale_fill_manual(name = "Mutation type", labels = c("missense", "nonsense", "synonymous"), values=c("#CC79A7", "#D55E00", "#0072B2")) +
  # geom_hline(yintercept=50, linetype="dashed", color = "black", size = 0.9) +
  ylab("7-mer substitution rate") +
  xlab("") +
  ggtitle("(A)") +
  ylim(0, 0.15) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.15)) +
  coord_flip() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = -2, size=6) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  # scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(colour = "black", size=1),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p_wilcox

### PLOT RATIOS

tmp <- unique(ct_psnv[,c("k7_mutation", "mutation", "Mutation_type")])
dt_counts <- as.data.table(table(tmp[,c("mutation", "Mutation_type")]))
dt_counts$fraction <- dt_counts$N / 4096
# dt_counts <- dt_counts[Mutation_type != "missense"]
dt_counts$Mutation_type <- factor(dt_counts$Mutation_type, levels = c("nonsense", "synonymous", "missense"))

p_fraction <- ggplot(data=dt_counts, aes(x=mutation, y=fraction, fill=Mutation_type)) +
  geom_bar(stat="identity", colour="black", position = "dodge") +
  scale_fill_manual(values=c("#D55E00", "#0072B2", "#CC79A7")) +
  coord_flip() +
  ggtitle("(C)") +
  xlab("") +
  ylab("Fraction of 7-mer substitutions with\npotential to cause mutation type") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p_fraction

### GRID ARRANGE

p1 <- grid.arrange(p_box, p_fraction, ncol = 2, nrow = 1)
p_out <- grid.arrange(p_wilcox, p1, ncol = 1, nrow = 2, heights = c(3,5))

ggsave(file = "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/figures_tables/figure_main_mutation_CDS_effect.jpg", p_out, height = 10, width = 8.5)


# dt_counts <- as.data.table(table(ct_psnv[,c("mutation", "Mutation_type")]))
# dt_counts <- dt_counts[Mutation_type != "missense"]
# dt_counts <- dt_counts[, total:=sum(N), by=Mutation_type]
# dt_counts$fraction <- dt_counts$N / dt_counts$total
# # dt_counts <- dt_counts[Mutation_type != "missense"]
#
# ggplot(data=dt_counts, aes(x=Mutation_type, y=fraction, fill=mutation)) +
#   geom_bar(stat="identity", colour="black") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         # legend.position = "top",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 14),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   )

#####

# p_box <- ggplot(ct_psnv, mapping = aes(x=mutation, y=k7_mu_rate)) +
#   # geom_violin(width=1) +
#   geom_boxplot(width=0.2, outlier.shape = NA, fill = "white") +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   # ggtitle("") +
#   ylim(0, 0.1) +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   )
# p_box
# 
# p_box_C <- ggplot(ct_psnv[mu_group == "C>"], mapping = aes(x=mutation, y=k7_mu_rate)) +
#   geom_violin(width=1.2) +
#   geom_boxplot(width=0.15, outlier.shape = NA, fill = "white") +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   ylim(0, 0.2) +
#   # ggtitle("") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# p_box_C
# 
# p_box_A <- ggplot(ct_psnv[mu_group == "A>"], mapping = aes(x=mutation, y=k7_mu_rate)) +
#   geom_violin(width=1.2) +
#   geom_boxplot(width=0.15, outlier.shape = NA, fill = "white") +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   # ggtitle("") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# p_box_A
# 
# p_box_CGTG <- ggplot(ct_psnv[mu_group == "CG>TG"], mapping = aes(x=mutation, y=k7_mu_rate)) +
#   # geom_violin(width=1.2) +
#   geom_boxplot(width=0.15, outlier.shape = NA, fill = "white") +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   # ggtitle("") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# p_box_CGTG
# 
# tmp <- ct_psnv[Mutation_type %in% c("nonsense", "synonymous", "missense") & mu_group %in% c("A>", "C>")]
# tmp$mu_group <- "(A|C)>"
# dt_plot <- rbind(tmp, ct_psnv[Mutation_type %in% c("nonsense", "synonymous", "missense") & mu_group %in% c("A>", "C>")])
# 
# # p_box <- ggplot(dt_plot, mapping = aes(x=mutation, y=k7_mu_rate, fill=Mutation_type)) +
# #   geom_boxplot(width=0.6, outlier.shape = NA, position=position_dodge(0.8)) +
# #   scale_fill_manual(values=c("#0072B2", "#D55E00")) +
# #   # geom_hline(yintercept=50, linetype="dashed", color = "black", size = 0.9) +
# #   ylab("7-mer substitution rate") +
# #   xlab("") +
# #   # ggtitle("(A) Distributions") +
# #   ylim(0, 0.1) +
# #   # scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.15)) +
# #   coord_flip() +
# #   theme_classic() +
# #   theme(axis.text = element_text(size = 14),
# #         axis.title = element_text(size = 14),
# #         plot.title =  element_text(size = 20, face = "bold"),
# #         # legend.position = "top",
# #         legend.title = element_blank(),
# #         legend.text = element_text(size = 14),
# #         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
# #         panel.border = element_rect(colour = "black", fill=NA, size=1)
# #   )
# # p_box
# 
# p_box <- ggplot(dt_plot, mapping = aes(x=mu_group, y=k7_mu_rate, fill=Mutation_type)) +
#   geom_boxplot(width=0.6, outlier.shape = NA, position=position_dodge(0.8)) +
#   # scale_fill_manual(values=c("#0072B2", "#D55E00")) +
#   # geom_hline(yintercept=50, linetype="dashed", color = "black", size = 0.9) +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   # ggtitle("(A) Distributions") +
#   ylim(0, 0.15) +
#   # scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.15)) +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         # legend.position = "top",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 14),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# p_box
# 
# dt_counts <- as.data.table(table(ct_psnv[,c("mutation", "Mutation_type")]))
# # dt_counts <- dt_counts[Mutation_type != "missense"]
# dt_counts <- dt_counts[, total:=sum(N), by=mutation]
# dt_counts$fraction <- dt_counts$N / dt_counts$total
# # dt_counts <- dt_counts[Mutation_type != "missense"]
# 
# ggplot(data=dt_counts, aes(x=mutation, y=N, fill=Mutation_type)) +
#   geom_bar(stat="identity", colour="black") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         # legend.position = "top",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 14),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# 
# # dt_counts <- as.data.table(table(ct_psnv[,c("mutation", "Mutation_type")]))
# # dt_counts <- dt_counts[Mutation_type != "missense"]
# # dt_counts <- dt_counts[, total:=sum(N), by=Mutation_type]
# # dt_counts$fraction <- dt_counts$N / dt_counts$total
# # # dt_counts <- dt_counts[Mutation_type != "missense"]
# # 
# # ggplot(data=dt_counts, aes(x=Mutation_type, y=fraction, fill=mutation)) +
# #   geom_bar(stat="identity", colour="black") +
# #   coord_flip() +
# #   theme_classic() +
# #   theme(axis.text = element_text(size = 14),
# #         axis.title = element_text(size = 14),
# #         plot.title =  element_text(size = 20, face = "bold"),
# #         # legend.position = "top",
# #         legend.title = element_blank(),
# #         legend.text = element_text(size = 14),
# #         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
# #         panel.border = element_rect(colour = "black", fill=NA, size=1)
# #   ) 