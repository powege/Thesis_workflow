rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(cowplot)

### SET VARS ==========

in.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Constraint_scores/WM_Harr_etal_2016_allSPECIES_constraint_scores.csv"

### IMPORT DATA ==========

dt <- fread(in.file)

### PLOT ==========

p_syn <- ggplot(dt, aes(x = exp_synonymous, y = n_synonymous)) +
  geom_point(aes(alpha = 0.1, colour = "grey")) +
  geom_abline(slope=1, intercept=0) +
  # scale_color_manual(breaks = c("Least constrained", "Most constrained"),
  #                    values=c('gray80', 'blue', 'red')) +
  ggtitle("(A)") +
  xlab("Expected synonymous variant sites") +
  ylab('Observed synonymous variant sites') +
  # xlim(0, 150) +
  # ylim(0, 150) +
  theme_bw() +
  theme(
    plot.title =  element_text(size = 20, face = "bold"),
    legend.position="none",
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    text = element_text(size = 14),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
p_syn

p_nsyn <- ggplot(dt, aes(x = exp_nonsynonymous, y = n_nonsynonymous)) +
  geom_point(aes(alpha = 0.1, colour = "grey")) +
  geom_abline(slope=1, intercept=0) +
  # scale_color_manual(breaks = c("Least constrained", "Most constrained"),
  #                    values=c('gray80', 'blue', 'red')) +
  ggtitle("(B)") +
  xlab("Expected nonsynonymous variant sites") +
  ylab('Observed nonsynonymous variant sites') +
  # xlim(0, 150) +
  # ylim(0, 150) +
  theme_bw() +
  theme(
    plot.title =  element_text(size = 20, face = "bold"),
    legend.position="none",
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    text = element_text(size = 14),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
p_nsyn

Fig_out <- plot_grid(p_syn, p_nsyn, ncol = 2, nrow = 1)

### EXPORT ==========
save_plot("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Figures_tables/figure_supplement_obs_exp.jpg", Fig_out, ncol = 1, nrow = 1, base_height = 5, base_width = 10)
