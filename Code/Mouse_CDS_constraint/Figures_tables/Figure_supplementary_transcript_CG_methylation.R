rm(list = ls())
graphics.off()

library(data.table)
library(cowplot)

### FUNCTIONS
source("~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R")

#### SET PATHS
pos.in.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_canPC_pos.csv"
meth.in.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz"
t.meth.in.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/m.mus_grc38_ensembl_transcript_methylation.csv"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Figures_tables/Figure_supplementary_transcript_CG_unmethylated.jpg"

### IMPORT 
dt <- fread(t.meth.in.file)
pos <- fread(pos.in.file)
meth_all <- fread(paste0("gunzip -cq ", meth.in.file))

### FORMAT

pos <- pos[ensembl_transcript_id %in% dt$ensembl_transcript_id]
pos <- pos[complete.cases(pos[,c("chromosome_name", "genomic_coding_start", "genomic_coding_end")])]

meth_all$chromosome <- as.character(meth_all$chromosome)

tran_CG <- bed.intersect(pos[,c("chromosome_name", "genomic_coding_start", "genomic_coding_end")],
                      meth_all[,c("chromosome", "start", "end")])

meth <- meth_all[tran_CG, on = colnames(tran_CG)]
meth$ES_coverage[is.na(meth$ES_coverage)] <- 0

dt_bar <- data.table(percent_methylated = meth$ES_percent_methylated[meth$ES_coverage >= 5])
dt_bar <- as.data.table(table(round(dt_bar$percent_methylated/5)*5))
dt_bar$percentage <- (dt_bar$N/sum(dt_bar$N))*100
dt_bar$V1 <- factor(dt_bar$V1, levels = dt_bar$V1)
dt_bar$group <- "unclassified"
dt_bar$group[dt_bar$V1 %in% seq(0, 20, by = 5)] <- "unmethylated"

ES_coverage <- nrow(meth[ES_coverage >= 5]) / nrow(meth)
percentage_unmethylated <- sum(dt_bar$percentage[dt_bar$group == "unmethylated"])
percentage_unclassified <- sum(dt_bar$percentage[dt_bar$group == "unclassified"])
percentage_unmethylated <- nrow(meth[ES_coverage >= 5 & ES_percent_methylated <= 20]) / nrow(meth)
percentage_unclassified <- 1 - (percentage_unmethylated)

dt$percentage_unmethylated <- ( dt$n_CGunmethylated / dt$n_CG )*100
dt_bar2 <- data.table(percent_unmethylated = dt$percentage_unmethylated)
dt_bar2 <- as.data.table(table(round(dt_bar2$percent_unmethylated/5)*5))
dt_bar2$percentage <- (dt_bar2$N/sum(dt_bar2$N))*100
dt_bar2$V1 <- factor(dt_bar2$V1, levels = dt_bar2$V1)

nrow(dt[n_CGunmethylated == 0]) / nrow(dt)


### PLOT 

fig_A <- ggplot(data=dt_bar, aes(x=V1, y=percentage, fill=group)) + 
  geom_bar(stat="identity", color="black") +
  xlab("percentage methylated") +
  ylab("percentage of CGs") +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  ggtitle("(A)") +
  theme_classic() +
  theme(
    plot.title =  element_text(size = 20, face = "bold"),
    legend.position = c(0.16, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.background = element_rect(linetype="solid", colour ="black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
fig_A

fig_B <- ggplot(data=dt_bar2, aes(x=V1, y=percentage)) + 
  geom_bar(stat="identity", color="black", fill = "#E69F00") +
  xlab("percentage of unmethylated CGs") +
  ylab("percentage of transcripts") +
  ggtitle("(B)") +
  theme_classic() +
  theme(
    plot.title =  element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
fig_B

out_plot <- plot_grid(fig_A, fig_B,
                      ncol = 1, nrow = 2)

### EXPORT
save_plot(filename = out.file,
          plot = out_plot,
          ncol = 1, nrow = 1,
          base_height = 8, base_width = 6)





                      
                      
                      
                      
                      
                      
                      
                      