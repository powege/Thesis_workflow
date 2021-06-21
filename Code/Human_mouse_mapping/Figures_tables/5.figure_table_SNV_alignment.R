### SCRIPT that creates tables and figures for ClinVar and GWAS SNV alignment
### INPUT:
# ClinVar alignment by chromosome
# GWAS alignnment by chromosome
### OUTPUT
# figure of GWAS SNV alignment by feature
# figure of ClinVar SNV alignment by feature
# table of GWAS and ClinVar SNV alignment by feature

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plyr)
library(gridExtra)
library(RColorBrewer)

### SET VARS
in.clinvar <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_SNV_alignment_by_chr.csv"
in.gwas <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_by_chr.csv"
in.gwas.null <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_NULL_by_chr.csv"
out.table <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Tables/Table_hsap_v_mmus_RegBuild_v101_SNV_alignment.csv"
out.clinvar.figure <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Figures/Figure_hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_SNV_alignment.jpeg"
out.gwas.figure <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/Figures/Figure_hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment.jpeg"

### IMPORT

cv_raw <- fread(in.clinvar)
gwas_raw <- fread(in.gwas)
gwas_null <- fread(in.gwas.null)

### FORMAT TABLE

# sub <- multi_raw[annotation == "Exon - CDS"]
grindilow <- function(sub){
  c(annotation = sub$annotation[1],
    apply(sub[,c("n_total", "n_aligned", "n_conserved")], 2, sum))
}
cv <- as.data.table(ddply(cv_raw, "annotation", grindilow))
gwas <- as.data.table(ddply(gwas_raw, "annotation", grindilow))

cv$n_conserved[cv$annotation %in% c("Total", "Functional - nonCDS", "Functional - all")] <- NA
gwas$n_conserved[gwas$annotation %in% c("Total", "Functional - nonCDS", "Functional - all")] <- NA

cv[, c(2:4)] <- lapply(cv[, c(2:4)], as.numeric)
gwas[, c(2:4)] <- lapply(gwas[, c(2:4)], as.numeric)

cv$percentage_aligned <- (cv$n_aligned / cv$n_total)*100
cv$percentage_conserved <- (cv$n_conserved / cv$n_total)*100
gwas$percentage_aligned <- (gwas$n_aligned / gwas$n_total)*100
gwas$percentage_conserved <- (gwas$n_conserved / gwas$n_total)*100

cv$group <- "ClinVar"
gwas$group <- "GWAS"
out_table <- rbind(cv, gwas)

### FORMAT NULL 

out_list_A <- list()
for(j in 1:length(unique(gwas_null$annotation))){
  sub <- gwas_null[annotation == unique(gwas_null$annotation)[j]]
  out_list_B <- list()
  for (i in 1:length(unique(sub$chromosome))){
    sub_sub <- sub[chromosome == unique(sub$chromosome)[i]]
    sub_sub$iteration <- 1:nrow(sub_sub)
    out_list_B[[i]] <- sub_sub
  }
  out_list_A[[j]] <- do.call("rbind", out_list_B)
}
gwas_null <- do.call("rbind", out_list_A)
gwas_null <- gwas_null[iteration <= 200]

out_list_A <- list()
for(j in 1:length(unique(gwas_null$annotation))){
  sub <- gwas_null[annotation == unique(gwas_null$annotation)[j]]
  out_list_B <- list()
  for (i in 1:length(unique(sub$iteration))){
    sub_sub <- sub[iteration == unique(sub$iteration)[i]]
    out_list_B[[i]] <- c(annotation = sub_sub$annotation[1],
                         iteration = sub_sub$iteration[1],
                         apply(sub_sub[,c("n_total", "n_aligned", "n_conserved")], 2, sum))
  }
  out_list_A[[j]] <- do.call("rbind", out_list_B)
}
gwas_null <- as.data.table(do.call("rbind", out_list_A))

gwas_null[, c(2:5)] <- lapply(gwas_null[, c(2:5)], as.numeric)
gwas_null$percentage_aligned_null <- (gwas_null$n_aligned / gwas_null$n_total)*100
gwas_null$percentage_conserved_null <- (gwas_null$n_conserved / gwas_null$n_total)*100

out_list <- list()
for(j in 1:length(unique(gwas_null$annotation))){
  sub <- gwas_null[annotation == unique(gwas_null$annotation)[j]]
  out_list[[j]] <- data.table(annotation = sub$annotation[i],
             mean_aligned_null = mean(sub$percentage_aligned_null),
             mean_conserved_null = mean(sub$percentage_conserved_null),
             sd_aligned_null = sd(sub$percentage_aligned_null),
             sd_conserved_null = sd(sub$percentage_conserved_null))
}
gwas_null <- do.call("rbind", out_list)
gwas_null$group <- "GWAS"

out_table <- gwas_null[out_table, on = c("annotation", "group")]

out_table$z_aligned <- (out_table$percentage_aligned - out_table$mean_aligned_null) / out_table$sd_aligned_null
out_table$z_conserved <- (out_table$percentage_conserved - out_table$mean_conserved_null) / out_table$sd_conserved_null
out_table$z_aligned_pval <- (1-pnorm(abs(out_table$z_aligned)))*2
out_table$z_conserved_pval <- (1-pnorm(abs(out_table$z_conserved)))*2

# for each iteration, plot the standard deviation






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
GW_cv_alignment <- cv$percentage_aligned[cv$annotation == "Total"]
GW_gwas_alignment <- gwas$percentage_aligned[gwas$annotation == "Total"]

cv_plot <- cv[!annotation %in% c("Functional - all", "Functional - nonCDS", "Total")]
cv_plot$x_lab <- paste0(cv_plot$annotation, "\n(n=", cv_plot$n_total, ")")
gwas_plot <- gwas[!annotation %in% c("Functional - all", "Functional - nonCDS", "Total")]
gwas_plot$x_lab <- paste0(gwas_plot$annotation, "\n(n=", gwas_plot$n_total, ")")

cv_plot$annotation <- factor(cv_plot$annotation, levels = as.character(order))
cv_plot <- cv_plot[order(annotation)]
lab_order <- cv_plot$x_lab
cv_plot$x_lab <- factor(cv_plot$x_lab, levels = as.character(lab_order))
gwas_plot$annotation <- factor(gwas_plot$annotation, levels = as.character(order))
gwas_plot <- gwas_plot[order(annotation)]
lab_order <- gwas_plot$x_lab
gwas_plot$x_lab <- factor(gwas_plot$x_lab, levels = as.character(lab_order))

cv_plot$percentage_not_conserved <- cv_plot$percentage_aligned - cv_plot$percentage_conserved
cv_plot <- melt(cv_plot, measure.vars = c("percentage_conserved", "percentage_not_conserved"))
cv_plot$variable <- factor(cv_plot$variable, levels = c("percentage_conserved", "percentage_not_conserved"))
cv_plot$col <- ifelse(cv_plot$variable == "percentage_not_conserved", "black", as.character(cv_plot$annotation))
cv_plot$col <- factor(cv_plot$col, levels = c("black", as.character(order)))
gwas_plot$percentage_not_conserved <- gwas_plot$percentage_aligned - gwas_plot$percentage_conserved
gwas_plot <- melt(gwas_plot, measure.vars = c("percentage_conserved", "percentage_not_conserved"))
gwas_plot$variable <- factor(gwas_plot$variable, levels = c("percentage_conserved", "percentage_not_conserved"))
gwas_plot$col <- ifelse(gwas_plot$variable == "percentage_not_conserved", "black", as.character(gwas_plot$annotation))
gwas_plot$col <- factor(gwas_plot$col, levels = c("black", as.character(order)))

cv_color_scale <- c("black",  "#B15928", brewer.pal(12, name = "Set3"))
cv_color_map <- setNames(cv_color_scale, levels(unique(cv_plot$col)))
gwas_color_scale <- c("black",  "#B15928", brewer.pal(12, name = "Set3"))
gwas_color_map <- setNames(gwas_color_scale, levels(unique(gwas_plot$col)))

### PLOT

p_cv <- ggplot(cv_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (% SNV sites)") +
  ggtitle("ClinVar Pathogenic") +
  scale_fill_manual(values = cv_color_map) +
  geom_hline(yintercept=GW_cv_alignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_cv

p_gwas <- ggplot(gwas_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (% SNV sites)") +
  ggtitle("GWAS DDS") +
  scale_fill_manual(values = gwas_color_map) +
  geom_hline(yintercept=GW_gwas_alignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_gwas

### EXPORT
ggsave(out.clinvar.figure, plot = p_cv, height = 6, width = 6)
ggsave(out.gwas.figure, plot = p_gwas, height = 6, width = 6)
fwrite(out_table, out.table)




