rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

# Correlation between regulatory feature constraint (UTRs, promoters, enhancers, CTCF binding sites) 
# and number of target genes.


### FUNCTION
# plot_file = "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_n_tss~annotation_constraint.csv"
# interest = c("CTCF binding","Enhancer - distal","Enhancer - proximal")
# y_min=0
# y_max=15
alakazam <- function(plot_file,
                     interest,
                     y_min,
                     y_max){
  
### IMPORT
dt_plot <- fread(plot_file)

### FORMAT

# subset annotations
dt_plot <- dt_plot[annotation %in% interest]

# aggregate by average tss
dt_plot_mn <- aggregate(list(mean_tss=dt_plot$n_tss), by=list(annotation=dt_plot$annotation, 
                                                              constraint_percentile=dt_plot$residual_percentile), 
                        FUN=mean)
dt_plot_md <- aggregate(list(mean_tss=dt_plot$n_tss), by=list(annotation=dt_plot$annotation, 
                                                              constraint_percentile=dt_plot$residual_percentile), 
                        FUN=median)

### BOXPLOT
# n_tss ~ annotation
box1 <- ggplot(dt_plot, aes(x=annotation, y=n_tss, fill=annotation)) +
  geom_boxplot() +
  coord_flip()
box1

### SCATTERPLOT
# n_tss ~ constraint percentile
p1 <- ggplot(dt_plot_mn, aes(x=constraint_percentile, y=mean_tss, colour=annotation)) +
  geom_point() +
  # geom_smooth(method = "lm") +
  geom_smooth(method = "loess") +
  ylab("Mean # target genes") +
  xlab("Constraint percentile rank") +
  # scale_y_continuous(breaks=c(1, 25, 50, 75, 100)) +
  # scale_x_continuous(breaks=c(1, 25, 50, 75, 100)) +
  ylim(y_min, y_max) +
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

return(p1)
}

### RUN 

m_plot <- alakazam(plot_file = "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_n_tss~annotation_constraint.csv",
                   interest = c("CTCF binding","Enhancer - distal","Enhancer - proximal"),
                   y_min=3,
                   y_max=12)
m_plot
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_M_n_tss~annotation_constraint.jpg", 
       plot = m_plot, height = 6, width = 8)

h_plot <- alakazam(plot_file = "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Human_n_tss~annotation_constraint.csv",
                   interest = c("CTCF binding","Enhancer - distal","Enhancer - proximal"),
                   y_min=3,
                   y_max=12)
h_plot
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_H_n_tss~annotation_constraint.jpg", 
       plot = h_plot, height = 6, width = 8)
#####

# Correlation in constraint between regulatory features (UTRs, promoters, enhancers, CTCF binding sites) 
# and target genes (oe ratio).

# Correlation in constraint between regulatory features (UTRs, promoters, enhancers, CTCF binding sites)
# of orthologous genes. 





