rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(ggplot2)

# Correlation between regulatory feature constraint (UTRs, promoters, enhancers, CTCF binding sites) 
# and nnumber of pathogenic SNV sites

### PLOT FOR MOUSE

### IMPORT 
dt_plot <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_SNVs_by_annotation_percentile.csv")

# aggregate by average OE
dt_plot_mn <- aggregate(list(mean_SNVs=dt_plot$fraction), by=list(annotation=dt_plot$annotation, 
                                                                constraint_percentile=dt_plot$percentile), 
                        FUN=mean)
dt_plot_mn$mean_SNVs <- dt_plot_mn$mean_SNVs*1000

### SCATTERPLOT
# n_SNVs ~ constraint percentile
pm <- ggplot(dt_plot_mn, aes(x=constraint_percentile, y=mean_SNVs, colour=annotation)) +
  geom_point(alpha = 1/3) +
  geom_smooth(method = "lm") +
  # geom_smooth(method = "loess") +
  ylab("Pathogenic SNV sites") +
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
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_n_SNV~M_annotation_constraint.jpg",
       plot = pm, height = 6, width = 8)


### PLOT FOR HUMAN



#####







