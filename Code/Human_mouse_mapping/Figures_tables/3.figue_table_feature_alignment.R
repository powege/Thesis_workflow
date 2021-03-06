### SCRIPT that creates tables and figures for feature alignment
### INPUT: 
# multicell feature alignment by chromosome
# tissue specific feature alignnment by chromosome
### OUTPUT
# table and figure of multicell feature alignment
# table of tissue specific feature alignment

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plyr)
library(gridExtra)
library(RColorBrewer)

### SET VARS 
in.multicell.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_summary.csv"
in.tissue.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_grch38_v_mmus_grcm38_v101_tissue_feature_alignment_summary.csv"
out.multicell.table <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Tables/Table_hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_summary.csv"
out.tissue.table <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Tables/Table_hsap_grch38_v_mmus_grcm38_v101_tissue_feature_alignment_summary.csv"
out.multicell.figure <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Figures/Figure_hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment.jpg"

### IMPORT
multi_raw <- fread(in.multicell.file)
tiss_raw <- fread(in.tissue.file)

### FORMAT TABLE 

# sub <- multi_raw[annotation == "Exon - CDS"]
grindilow <- function(sub){
  c(annotation = sub$annotation[1], 
    apply(sub[,c("n_total", "n_aligned", "n_conserved")], 2, sum))
}
multi <- as.data.table(ddply(multi_raw, "annotation", grindilow))

grindilow2 <- function(sub){
  c(tissue = sub$tissue[1], 
    annotation = sub$annotation[1], 
    apply(sub[,c("n_total", "n_aligned", "n_conserved")], 2, sum))
}
tiss <- as.data.table(rbind(ddply(tiss_raw[tissue == "heart"], "annotation", grindilow2),
               ddply(tiss_raw[tissue == "kidney"], "annotation", grindilow2),
               ddply(tiss_raw[tissue == "spleen"], "annotation", grindilow2)))

multi$n_conserved[multi$annotation %in% c("Total", "Functional - nonCDS", "Functional - all")] <- NA

multi[, c(2:4)] <- lapply(multi[, c(2:4)], as.numeric)
tiss[, c(3:5)] <- lapply(tiss[, c(3:5)], as.numeric)

multi$percentage_aligned <- round( (multi$n_aligned / multi$n_total)*100, 2)
multi$percentage_conserved <- round( (multi$n_conserved / multi$n_total)*100, 2)
tiss$percentage_aligned <- round( (tiss$n_aligned / tiss$n_total)*100, 2)
tiss$percentage_conserved <- round( (tiss$n_conserved / tiss$n_total)*100, 2)
multi$genomic_coverage <- round( (multi$n_total/multi$n_total[multi$annotation == "Total"])*100, 2)
tiss$genomic_coverage <- round( (tiss$n_total/multi$n_total[multi$annotation == "Total"])*100, 2)

### FORMAT FIGURES

order <- c("Exon - CDS",
           "Exon - 5'UTR",
           "Exon - 3'UTR",
           "Exon - other",
           "Promoter",
           "Intron - proximal",
           "Enhancer - proximal",
           "Enhancer - distal",
           "CTCF binding",
           "TAD boundry",
           "Miscellaneous",
           "Intron - distal",
           "Unannotated")
order <- rev(order)
GW_Halignment <- multi$percentage_aligned[multi$annotation == "Total"]

multi_plot <- multi[!annotation %in% c("Functional - all", "Functional - nonCDS", "Total")]
multi_plot$x_lab <- paste0(multi_plot$annotation, "\n(", multi_plot$genomic_coverage, ")")

multi_plot$annotation <- factor(multi_plot$annotation, levels = as.character(order))
multi_plot <- multi_plot[order(annotation)]
lab_order <- multi_plot$x_lab
multi_plot$x_lab <- factor(multi_plot$x_lab, levels = as.character(lab_order))

multi_plot$percentage_not_conserved <- multi_plot$percentage_aligned - multi_plot$percentage_conserved
multi_plot <- melt(multi_plot, measure.vars = c("percentage_conserved", "percentage_not_conserved"))
multi_plot$variable <- factor(multi_plot$variable, levels = c("percentage_conserved", "percentage_not_conserved"))
multi_plot$col <- ifelse(multi_plot$variable == "percentage_not_conserved", "black", as.character(multi_plot$annotation))
multi_plot$col <- factor(multi_plot$col, levels = c("black", as.character(order)))

color_scale <- c("black",  "#B15928", brewer.pal(12, name = "Set3"))
color_map <- setNames(color_scale, levels(unique(multi_plot$col)))

### PLOT

p_multi <- ggplot(multi_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (% bp)") +
  # ggtitle("B") +
  scale_fill_manual(values = color_map) +
  geom_hline(yintercept=GW_Halignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_multi

### EXPORT
ggsave(out.multicell.figure, plot = p_multi, height = 6, width = 6)
fwrite(multi, out.multicell.table)
fwrite(tiss, out.tissue.table)


