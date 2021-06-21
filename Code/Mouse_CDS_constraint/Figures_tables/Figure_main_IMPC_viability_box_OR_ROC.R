rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plotROC)

### PLOT OR 

### IMPORT 
dt_or <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_constraint_score_IMPC_viability_fertility_OR.csv")

### FORMAT
scores <- c("OE_nonsynonymous_percentile", "spretus_dNdS_percentile")
groups <- c("L", "SV", "VP", "VN")

dt_or <- dt_or[score %in% scores & group %in% groups]
dt_or$CI_min <- dt_or$OR-dt_or$CI95
dt_or$CI_max <- dt_or$OR+dt_or$CI95
dt_or$group <- ordered(dt_or$group, levels = c("L", "SV", "VP", "VN"))
dt_or$score[dt_or$score == "OE_nonsynonymous_percentile"] <- "NOER"
dt_or$score[dt_or$score == "spretus_dNdS_percentile"] <- "dNdS"
dt_or_high <- dt_or[threshold == ">=91"]
dt_or_low <- dt_or[threshold == "<=10"]

# plot most constraied OR
p_or_high <- ggplot(dt_or_high, aes(x=group, y=OR, color=score)) + 
  geom_errorbar(aes(ymax = CI_max, ymin = CI_min), 
                stat = "identity",
                position = position_dodge(width = 0.5),
                width=0.2,
                size=1.2) + 
  geom_point(size = 2.5,
             position = position_dodge(width = 0.5)) + 
  scale_color_manual(values=c("#0072B2", "#D55E00")) +
  xlab("") +
  ylab("Odds ratio") +
  ggtitle("(C) 10% most constrained") +
  # scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  #                    limits = c(0, 1)) +
  coord_flip() +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size = 0.9) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p_or_high

# plot least constrained OR
p_or_low <- ggplot(dt_or_low, aes(x=group, y=OR, color=score)) + 
  geom_errorbar(aes(ymax = CI_max, ymin = CI_min), 
                stat = "identity",
                position = position_dodge(width = 0.5),
                width=0.2,
                size=1.2) + 
  geom_point(size = 2.5,
             position = position_dodge(width = 0.5)) + 
  scale_color_manual(values=c("#0072B2", "#D55E00")) +
  xlab("") +
  ylab("Odds ratio") +
  ggtitle("(D) 10% least constrained") +
  # scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  #                    limits = c(0, 1)) +
  coord_flip() +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size = 0.9) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        legend.position = c(0.8, 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p_or_low


### PLOT distributions

### IMPORT 
dt_box <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv")

### FORMAT
dt_box <- dt_box[,c("mmus_external_gene_name",
                      "OE_nonsynonymous",
                      "spretus_dNdS",
                      "viability_status")]
# dt_box <- dt_box[complete.cases(dt_box) & viability_status != "",]
dt_box <- dt_box[viability_status != "",]

# percentiles
percentile <- ecdf(dt_box$OE_nonsynonymous[!duplicated(dt_box$mmus_external_gene_name)])
dt_box$OE_nonsynonymous_percentile <- 101 - (ceiling(percentile(dt_box$OE_nonsynonymous)*100))
percentile <- ecdf(dt_box$spretus_dNdS[!duplicated(dt_box$mmus_external_gene_name)])
dt_box$spretus_dNdS_percentile <- 101 - (ceiling(percentile(dt_box$spretus_dNdS)*100))

dt_long <- melt(dt_box,
                  id.vars=c("mmus_external_gene_name", "viability_status"), # ID variables - all the variables to keep but not split apart on
                  measure.vars=c("OE_nonsynonymous_percentile", "spretus_dNdS_percentile") # The source columns
                  # variable.name="condition", # Name of the destination column that will identify the original
                  # value.name="measurement" # column that the measurement came from
)
dt_long$viability_status <- ordered(dt_long$viability_status, levels = c("L", "SV", "VP", "VN"))
dt_long$variable <- as.character(dt_long$variable)
dt_long$variable[dt_long$variable == "OE_nonsynonymous_percentile"] <- "NOER"
dt_long$variable[dt_long$variable == "spretus_dNdS_percentile"] <- "dNdS"

# plot 
p_box <- ggplot(dt_long, mapping = aes(x=viability_status, y=value, group=interaction(viability_status, variable))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = 90, ymax = 100), fill = "grey", color = NA, alpha = 0.5, show.legend = F) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = 1, ymax = 10), fill = "grey", color = NA, alpha = 0.5, show.legend = F) +
  geom_violin(width=1, position=position_dodge(0.75), bw=1.5, aes(fill = variable)) +
  geom_boxplot(width=0.2, outlier.shape = NA, position=position_dodge(0.75), fill = "white") +
  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
  geom_hline(yintercept=50, linetype="dashed", color = "black", size = 0.9) +
  ylab("Constraint percentile rank") +
  xlab("") +
  ggtitle("(A) Distributions") +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100)) +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        # legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) 
p_box


### PLOT ROC 

### IMPORT 
dt_roc_raw <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_constraint_score_IMPC_viability_prediction.csv")
dt_roc <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_constraint_score_IMPC_viability_fertility_ROC.csv")

### FORMAT
dt_roc_raw <- dt_roc_raw[,c("mmus_external_gene_name", "obs", "OE_nonsynonymous_percentile_pred", "spretus_dNdS_percentile_pred")]

auc.OE <- round(dt_roc$ROC[dt_roc$group == "L_SV" & dt_roc$variable == "OE_nonsynonymous_percentile"], 2)
auc.dNdS <- round(dt_roc$ROC[dt_roc$group == "L_SV" & dt_roc$variable == "spretus_dNdS_percentile"], 2)
p.OE <- format(dt_roc$p_val[dt_roc$group == "L_SV" & dt_roc$variable == "OE_nonsynonymous_percentile"], digits = 2)
p.dNdS <- format(dt_roc$p_val[dt_roc$group == "L_SV" & dt_roc$variable == "spretus_dNdS_percentile"], digits = 2)

colnames(dt_roc_raw) <- c("mmus_external_gene_name", "obs",
                       paste0("NOER\nAUC=", auc.OE, " p=", p.OE),
                       paste0("dNdS\nAUC=", auc.dNdS, " p=", p.dNdS))

dt_roc_raw <- melt(dt_roc_raw, id.vars = c("mmus_external_gene_name", "obs"))
dt_roc_raw$variable <- as.factor(dt_roc_raw$variable)


### PLOT
p_roc <- ggplot(dt_roc_raw, aes(d = obs, m = value, color = variable)) + 
  geom_roc(n.cuts = 0) +
  style_roc() +
  scale_colour_manual(values=c("#D55E00", "#0072B2")) +
  ggtitle("(B) ROC (L and SV)") +
  # theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        legend.position = c(0.7, 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) 
  p_roc

  
  ### OUTPUT
  
  pout <- grid.arrange(p_box, p_roc, p_or_high, p_or_low, nrow = 2, ncol = 2, widths = c(5, 5))
  ggsave("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Figures_tables/figure_main_IMPC_viability_box_OR_ROC.jpg", plot = pout, height = 11, width = 11)
  


#####

  
  
  
  
