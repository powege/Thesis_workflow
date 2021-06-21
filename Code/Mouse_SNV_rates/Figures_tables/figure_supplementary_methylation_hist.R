rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)

### IMPORT 
meth <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz")

### FORMAT

ES_coverage <- nrow(meth[ES_coverage != 0]) / nrow(meth)
ES_coverage <- nrow(meth[ES_coverage >= 5]) / nrow(meth)

dt_bar <- data.table(percent_methylated = meth$ES_percent_methylated[meth$ES_coverage >= 5])
dt_bar <- as.data.table(table(round(dt_bar$percent_methylated/5)*5))
dt_bar$percentage <- (dt_bar$N/sum(dt_bar$N))*100
dt_bar$V1 <- factor(dt_bar$V1, levels = dt_bar$V1)
dt_bar$group <- "unclassified"
dt_bar$group[dt_bar$V1 %in% seq(0, 20, by = 5)] <- "unmethylated"
dt_bar$group[dt_bar$V1 %in% seq(60, 100, by = 5)] <- "methylated"

fig_bar <- ggplot(data=dt_bar, aes(x=V1, y=percentage, fill=group)) + 
  geom_bar(stat="identity", color="black") +
  xlab("percentage methylated") +
  ylab("percentage of CGs") +
  scale_fill_manual(values = c("#56B4E9", "#999999", "#E69F00")) +
  theme_classic() +
  theme(
    legend.position = c(0.16, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.background = element_rect(linetype="solid", colour ="black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
fig_bar

percentage_methylated <- sum(dt_bar$percentage[dt_bar$group == "methylated"])
percentage_unmethylated <- sum(dt_bar$percentage[dt_bar$group == "unmethylated"])
percentage_unclassified <- sum(dt_bar$percentage[dt_bar$group == "unclassified"])

percentage_methylated <- nrow(meth[ES_coverage >= 5 & ES_percent_methylated >= 60]) / nrow(meth)
percentage_unmethylated <- nrow(meth[ES_coverage >= 5 & ES_percent_methylated <= 20]) / nrow(meth)
percentage_unclassified <- 1 - (percentage_methylated + percentage_unmethylated)


### Export ==========
ggsave("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/figures_tables/figure_supplementary_methylation_hist.jpg", 
       fig_bar, 
       width = 6.5, 
       height = 3.5)



#####

# fig_bar <- ggplot(data=dt_bar, aes(x=V1, y=percentage)) + 
#   geom_bar(stat="identity", color="black", fill="grey") +
#   geom_rect(aes(xmin = 0.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
#             fill = "transparent", color = "red", size = 1) +
#   geom_rect(aes(xmin = 12.5, xmax = 21.5, ymin = -Inf, ymax = Inf),
#             fill = "transparent", color = "blue", size = 1) +
#   xlab("percentage methylated") +
#   ylab("percentage of CGs") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   )
# fig_bar


  

 





