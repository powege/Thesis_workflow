# Label Brca1 and Brca, Msh1 and Mlh2

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(scales)

### FUNCTIONS 

### FUNCTION for plotting Spearman's correlation on ggplot
# k <- cor.test(df_md_all$x, df_md_all$Category, method = "spearman")
rs_corr_eqn <- function(k, digits) {
  output <- substitute(italic(rho)~"="~rhovalue*","~italic(p)~"="~pvalue, 
                       list(rhovalue = unname(format(k$estimate, digits = digits)), 
                            pvalue = unname(format(k$p.value, digits = digits))))
  as.character(as.expression(output))  
}

### SET VARS
infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv"

### IMPORT
dt <- fread(infile)

### PLOT COR MATRIX
dt_plot <- cor(dt[,c("OE_nonsynonymous_percentile",
                     # "Z_nonsynonymous_percentile",
                     "spretus_dNdS_percentile",
                     "mis_z_percentile",
                     "oe_lof_upper_percentile",
                     "pLI_percentile")], 
               method = "spearman", use = "na.or.complete")
dt_plot <- as.matrix(dt_plot)
colnames(dt_plot) <- c(  "Mouse OE",
                         # "Mouse nonsynonymous z-score",
                         "Mouse dN/dS",
                         "Missense z-score",
                         "LOEUF",
                         "pLI")
rownames(dt_plot) <- c(  "Mouse OE",
                         # "Mouse nonsynonymous z-score",
                         "Mouse dN/dS",
                         "Missense z-score",
                         "LOEUF",
                         "pLI")
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(dt_plot)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat$value <- round(melted_cormat$value, digits = 2)
melted_cormat$value <- abs(melted_cormat$value)

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "#D55E00", high = "#56B4E9", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Rank\nCorrelation") +
  ggtitle("(A)") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1),
        axis.text.y = element_text(vjust = 1, 
                                   size = 14, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value, fontface=2), color = "black", size = 5) +
  theme(
    plot.title =  element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(1.3, 0.0),
    # legend.direction = "vertical",
    plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
    text = element_text(size = 14)
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 12,
                               title.position = "top", title.hjust = 0.5))
ggheatmap 
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

### PLOT CLINVAR 

### FORMAT 

dt_clinvar <- dt[!is.na(dt$cds_length)]
dt_clinvar <- dt_clinvar[hsap_orthology_type == "ortholog_one2one"]
dt_clinvar <- dt_clinvar[!which(is.na(dt_clinvar$n_pathogenic) & is.na(dt_clinvar$n_benign)),]
dt_clinvar$n_pathogenic[is.na(dt_clinvar$n_pathogenic)] <- 0
dt_clinvar$n_benign[is.na(dt_clinvar$n_benign)] <- 0
dt_clinvar$n_pathogenic_kb <- (dt_clinvar$n_pathogenic / dt_clinvar$cds_length) * 1000 
dt_clinvar$n_benign_kb <- (dt_clinvar$n_benign / dt_clinvar$cds_length) * 1000 
dt_clinvar <- dt_clinvar[,c("mmus_external_gene_name",
               "OE_nonsynonymous",
               "OE_nonsynonymous_percentile",
               "n_pathogenic_kb",
               "n_benign_kb"
               )]
dt_clinvar <- dt_clinvar[complete.cases(dt_clinvar)]
# percentile <- ecdf(dt_clinvar$OE_nonsynonymous[!duplicated(dt_clinvar$mmus_external_gene_name)])
# dt_clinvar$OE_nonsynonymous_percentile <- 101 - (ceiling(percentile(dt_clinvar$OE_nonsynonymous)*100))

# hist(dt_clinvar$OE_nonsynonymous_percentile)
dt_pathogenic <- aggregate(list(average_pathogenic=dt_clinvar$n_pathogenic_kb,
                                average_benign=dt_clinvar$n_benign_kb), 
                           by=list(OE_nonsynonymous_percentile=dt_clinvar$OE_nonsynonymous_percentile), 
                           FUN=mean)

pathogenic_cor <- rs_corr_eqn(cor.test(dt_pathogenic$OE_nonsynonymous_percentile, dt_pathogenic$average_pathogenic), 2)
benign_cor <- rs_corr_eqn(cor.test(dt_pathogenic$OE_nonsynonymous_percentile, dt_pathogenic$average_benign), 2)

dt_pathogenic <- melt(dt_pathogenic, id.vars = "OE_nonsynonymous_percentile", easure.vars = c("average_pathogenic", "average_benign"))
dt_pathogenic$variable <- as.character(dt_pathogenic$variable)
dt_pathogenic$variable[dt_pathogenic$variable == "average_pathogenic"] <- "Pathogenic\nrho=0.27, p=0.0075"
dt_pathogenic$variable[dt_pathogenic$variable == "average_benign"] <- "Benign\nrho=-0.59, p=2.4e-10"

p_clinvar <- ggplot(dt_pathogenic, aes(x = OE_nonsynonymous_percentile, y = value, colour = variable)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method=lm, se=F) +
  scale_colour_manual(values=c("#D55E00", "#0072B2")) +
  # annotate("text", x = 37, y = max(dt_pathogenic$value), label = pathogenic_cor, colour="#D55E00", size = 5, parse=TRUE) +
  # annotate("text", x = 77, y = max(dt_pathogenic$value), label = benign_cor, colour="#0072B2", size = 5, parse=TRUE) +
  xlab("OE percentile rank") +
  ylab('mean SNVs per kb\nin the human orthologue') +
  ggtitle("(B)") +
  theme_bw() +
  theme(
        plot.title =  element_text(size = 20, face = "bold"),
        text = element_text(size = 14),
        # plot.title =  element_text(size = 20, face = "bold"),
        legend.position = c(0.65, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))
p_clinvar


### OUTPUT

pout <- grid.arrange(ggheatmap, p_clinvar, nrow = 1, ncol = 2, widths = c(5.5, 4.5))
ggsave("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Figures_tables/figure_main_score_matrix_Clinvar.jpg", plot = pout, height = 5, width = 10)


#####

# x <- dt[,c("OE_nonsynonymous_percentile",
#            "spretus_dNdS_percentile",
#            "mis_z_percentile",
#            "oe_lof_upper_percentile",
#            "pLI_percentile")]
# x <- x[complete.cases(x),]
# y <- cor.test(x$OE_nonsynonymous_percentile, x$spretus_dNdS_percentile,
#                method = "spearman", use = "na.or.complete")

# sum(dt_clinvar$n_pathogenic)
# sum(dt_clinvar$n_benign)
# length(unique(dt_clinvar$hsap_ensembl_transcript_id))



