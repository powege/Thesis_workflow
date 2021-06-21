rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(ggplot2)

# Correlation in constraint between regulatory features (UTRs, promoters, enhancers, CTCF binding sites) 
# and closest gene (oe ratio).
# Plot fold change relative to least constrained percentile

### PLOT FOR MOUSE

### IMPORT 
dt_plot <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_gene_constraint~annotation_constraint.csv")

# aggregate by average OE
dt_plot_mn <- aggregate(list(mean_OE=dt_plot$H_oe_lof_upper_percentile), by=list(annotation=dt_plot$M_annotation, 
                                                                                 constraint_percentile=dt_plot$M_residual_percentile), 
                        FUN=mean)
dt_plot_md <- aggregate(list(mean_OE=dt_plot$H_oe_lof_upper_percentile), by=list(M_annotation=dt_plot$M_annotation, 
                                                                                 M_constraint_percentile=dt_plot$M_residual_percentile), 
                        FUN=median)


### SCATTERPLOT
# human OE ~ constraint percentile
pm <- ggplot(dt_plot_mn, aes(x=constraint_percentile, y=mean_OE, colour=annotation)) +
  geom_point(alpha = 1/3) +
  geom_smooth(method = "lm") +
  # geom_smooth(method = "loess") +
  ylab("Closest genes constraint") +
  xlab("Constraint percentile rank") +
  # scale_y_continuous(breaks=c(1, 25, 50, 75, 100)) +
  # scale_x_continuous(breaks=c(1, 25, 50, 75, 100)) +
  # ylim(y_min, y_max) +
  xlim(0, 100) +
  scale_fill_brewer(palette="Set2") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        plot.margin=unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 
pm
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_gene~constraint~M_annotation_constraint.jpg", 
       plot = pm, height = 6, width = 8)


### PLOT FOR HUMAN

### IMPORT 
dt_plot <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Human_gene_constraint~annotation_constraint.csv")

# aggregate by average OE
dt_plot_mn <- aggregate(list(mean_OE=dt_plot$H_oe_lof_upper_percentile), by=list(annotation=dt_plot$H_annotation, 
                                                                                 constraint_percentile=dt_plot$H_residual_percentile), 
                        FUN=mean)
dt_plot_md <- aggregate(list(mean_OE=dt_plot$H_oe_lof_upper_percentile), by=list(M_annotation=dt_plot$H_annotation, 
                                                                                 M_constraint_percentile=dt_plot$H_residual_percentile), 
                        FUN=median)


### SCATTERPLOT
# human OE ~ constraint percentile
ph <- ggplot(dt_plot_mn, aes(x=constraint_percentile, y=mean_OE, colour=annotation)) +
  geom_point(alpha = 1/3) +
  geom_smooth(method = "lm") +
  # geom_smooth(method = "loess") +
  ylab("Closest genes constraint") +
  xlab("Constraint percentile rank") +
  # scale_y_continuous(breaks=c(1, 25, 50, 75, 100)) +
  # scale_x_continuous(breaks=c(1, 25, 50, 75, 100)) +
  # ylim(y_min, y_max) +
  xlim(0, 100) +
  scale_fill_brewer(palette="Set2") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        plot.margin=unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ph
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_gene~constraint~H_annotation_constraint.jpg", 
       plot = ph, height = 6, width = 8)

#####







