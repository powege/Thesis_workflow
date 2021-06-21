rm(list=ls())
graphics.off()

library(data.table)
library(plyr)
library(Hmisc)

### IMPORT
dt <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv")

### FORMAT

dt_plot <- cor(dt[,c("OE_nonsynonymous_percentile",
                      "Z_nonsynonymous_percentile",
                      "spretus_dNdS_percentile",
                      "mis_z_percentile",
                      "oe_lof_upper_percentile",
                      "pLI_percentile")], 
               method = "spearman", use = "na.or.complete")
dt_plot <- as.matrix(dt_plot)
colnames(dt_plot) <- c(  "Mouse nonsynonymous OE",
                         "Mouse nonsynonymous z-score",
                         "Mouse dN/dS",
                         "Human missense z-score",
                         "Human LOEUF",
                         "Human pLI")
rownames(dt_plot) <- c(  "Mouse nonsynonymous OE",
                         "Mouse nonsynonymous z-score",
                         "Mouse dN/dS",
                         "Human missense z-score",
                         "Human LOEUF",
                         "Human pLI")
melted_cormat <- melt(dt_plot, na.rm = TRUE)
melted_cormat$value <- round(melted_cormat$value, digits = 2)

### EXPORT

# fwrite(melted_cormat, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/constraint_score_correlation.csv")

