rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

# Correlation in constraint between regulatory features (UTRs, promoters, enhancers, CTCF binding sites)
# of orthologous genes. 


### IMPORT 
dt_plot <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/H_annotation_constraint~M_annotation_constraint.csv")

# aggregate by average OE
dt_plot_mn <- aggregate(list(mean_H_constraint=dt_plot$H_residual_percentile), by=list(M_annotation=dt_plot$M_annotation, 
                                                                                 M_constraint_percentile=dt_plot$M_residual_percentile), 
                        FUN=mean)
dt_plot_md <- aggregate(list(mean_H_constraint=dt_plot$H_residual_percentile), by=list(M_annotation=dt_plot$M_annotation, 
                                                                                 M_constraint_percentile=dt_plot$M_residual_percentile), 
                        FUN=median)


### SCATTERPLOT
# human OE ~ constraint percentile
p1 <- ggplot(dt_plot_mn, aes(x=M_constraint_percentile, y=mean_H_constraint, colour=M_annotation)) +
  geom_point(alpha = 1/3) +
  geom_smooth(method = "lm") +
  # geom_smooth(method = "loess") +
  ylab("Human constraint") +
  xlab("Mouse constraint") +
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
p1
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_H_annotation~constraint~M_annotation_constraint.jpg", 
       plot = p1, height = 6, width = 8)


#####






