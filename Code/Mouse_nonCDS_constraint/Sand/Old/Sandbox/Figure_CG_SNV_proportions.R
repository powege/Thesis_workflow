# SCRIPT that plots:
#   the proportion of CG dinucleotides by annotation
#   the proportion of SNVs by annotation 

# change to bar plot with totals !!!!

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)

### FUNCTIONS

alakazam <- function(hvar_file,
                     mvar_file,
                     h_av,
                     m_av,
                     order){

### IMPORT
hvar <- fread(hvar_file)
mvar <- fread(mvar_file)

### FORMAT

# colnames
colnames(hvar) <- c("chromosome", "start", "end", "annotation", "strand", "ID", "n_var")
colnames(mvar) <- c("chromosome", "start", "end", "annotation", "strand", "ID", "n_var")

# calculate proportions by sequence
hvar$p_var <- hvar$n_var / ((hvar$end +1) - hvar$start)
mvar$p_var <- mvar$n_var / ((mvar$end +1) - mvar$start)

# # calculate average
# h_av <- median(hvar$p_var)
# m_av <- median(mvar$p_var)

# calculate deviation from average
hvar$p_var_dif <- hvar$p_var - h_av
mvar$p_var_dif <- mvar$p_var - m_av

# add species
hvar$species <- "Human"
mvar$species <- "Mouse"

# rbind
dt <- rbind(hvar, mvar)
rm(hvar, mvar)

# subset annotations
dt <- dt[annotation %in% order]

# calculate 1st 2nd and 3rd quartiles for proportions by annotation
tmp1 <- aggregate(dt$p_var, list(dt$species, dt$annotation), median)
colnames(tmp1) <- c("Species", "Category", "Q2")
tmp2 <- aggregate(dt$p_var, list(dt$species, dt$annotation), quantile, probs = 0.25)
colnames(tmp2) <- c("Species", "Category", "Q1")
tmp3 <- aggregate(dt$p_var, list(dt$species, dt$annotation), quantile, probs = 0.75)
colnames(tmp3) <- c("Species", "Category", "Q3")
q_prop <- merge(tmp2, tmp1)
q_prop <- merge(q_prop, tmp3)
rm(tmp1, tmp2, tmp3)

# calculate fold change B/A - 1
q_fc <- q_prop
q_fc$Q1[q_fc$Species == "Mouse"] <- (q_fc$Q1[q_fc$Species == "Mouse"]/m_av)
q_fc$Q2[q_fc$Species == "Mouse"] <- (q_fc$Q2[q_fc$Species == "Mouse"]/m_av)
q_fc$Q3[q_fc$Species == "Mouse"] <- (q_fc$Q3[q_fc$Species == "Mouse"]/m_av)
q_fc$Q1[q_fc$Species == "Human"] <- (q_fc$Q1[q_fc$Species == "Human"]/h_av)
q_fc$Q2[q_fc$Species == "Human"] <- (q_fc$Q2[q_fc$Species == "Human"]/h_av)
q_fc$Q3[q_fc$Species == "Human"] <- (q_fc$Q3[q_fc$Species == "Human"]/h_av)

# set rectangle coordinates 
rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1), 
                    xmax = tail(seq, -1), 
                    Category = levels(as.factor(order)),
                    rect_type = c("a", "c"))
rects <- rbind(rects)
q_prop <- merge(q_prop, rects)
q_fc <- merge(q_fc, rects)

# set factor order
q_prop$Category <- factor(q_prop$Category, levels = order)
q_fc$Category <- factor(q_prop$Category, levels = order)

# output
output <- list(dt, q_prop, q_fc)
return(output)
}

### IMPORT

# CG
h_cg_file <- "~/Dropbox/PhD/Data/NC_constraint/Variables/Human_GRC38_GENCODE_RegBuild_annotation_nCG.csv"
m_cg_file <- "~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_GRC38_GENCODE_RegBuild_annotation_nCG.csv"

# SNV
h_snv_file <- "~/Dropbox/PhD/Data/NC_constraint/Variables/Human_GRC38_GENCODE_RegBuild_annotation_nSNV.csv"
m_snv_file <- "~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_GRC38_GENCODE_RegBuild_annotation_nSNV.csv"

### SET VARS
order <-  rev(c("Exon - CDS", 
            "Exon - UTR",
            "Exon - other",
            "Promoter", 
            "Enhancer - proximal", 
            "Enhancer - distal", 
            "CTCF binding", 
            "Miscellaneous", 
            "Intron - distal",
            "Unannotated"))

### FORMAT

# CG
plot_cg_all <- alakazam(hvar_file = h_cg_file, 
                        mvar_file = m_cg_file, 
                        h_av = 0.052,
                        m_av = 0.052,
                        order)
plot_cg <- plot_cg_all[[1]]
plot_cg_prop <- plot_cg_all[[2]]
plot_cg_prop$Q1 <- plot_cg_prop$Q1 * 1000
plot_cg_prop$Q2 <- plot_cg_prop$Q2 * 1000
plot_cg_prop$Q3 <- plot_cg_prop$Q3 * 1000
plot_cg_fc <- plot_cg_all[[3]]

# SNV
plot_snv_all <- alakazam(hvar_file = h_snv_file, 
                        mvar_file = m_snv_file, 
                        h_av = 0.007,
                        m_av = 0.012,
                        order)
plot_snv <- plot_snv_all[[1]]
plot_snv_prop <- plot_snv_all[[2]]
plot_snv_prop$Q1 <- plot_snv_prop$Q1 * 1000
plot_snv_prop$Q2 <- plot_snv_prop$Q2 * 1000
plot_snv_prop$Q3 <- plot_snv_prop$Q3 * 1000
plot_snv_fc <- plot_snv_all[[3]]

### PLOT

# CG
pCG_prop <- ggplot(plot_cg_prop, aes(x=Category, y=Q2, color=Species)) + 
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.5,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  geom_rect(
    aes(xmin = plot_cg_prop$xmin,
        xmax = plot_cg_prop$xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = plot_cg_prop$rect_type),
    color = NA,
    alpha = 0.5,
    show.legend = F) +
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.4,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  # geom_hline(yintercept=1,
  #            linetype="dashed",
  #            color = "black",
  #            size=1) +
  # geom_hline(yintercept=m_av_SNV,
  #            linetype="dashed",
  #            color = "blue",
  #            size=1) +
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Genomic annotation") +
  ylab("CG dinucleotides per kb") +
  # ggtitle("CG dinucleotides") +
  # scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2),
  #                    limits = c(0, 2)) +
  # scale_x_discrete(position = "top") +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "none",
    # legend.title = element_blank(),
    # legend.justification=c(1,0),
    # legend.position=c(0.23, 0.70),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  )
pCG_prop

# pout <- grid.arrange(p1, p2, nrow = 1, widths = c(6, 6))
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_ann_CpG_SNV_prop.jpg", plot = pout, height = 6, width = 10)

# SNV
pSNV_prop <- ggplot(plot_snv_prop, aes(x=Category, y=Q2, color=Species)) + 
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.5,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  geom_rect(
    aes(xmin = plot_snv_prop$xmin,
        xmax = plot_snv_prop$xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = plot_snv_prop$rect_type),
    color = NA,
    alpha = 0.5,
    show.legend = F) +
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.4,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  # geom_hline(yintercept=1,
  #            linetype="dashed",
  #            color = "black",
  #            size=1) +
  # geom_hline(yintercept=m_av_SNV,
  #            linetype="dashed",
  #            color = "blue",
  #            size=1) +
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Genomic annotation") +
  ylab("SNV sites per kb") +
  # ggtitle("SNV sites") +
  # scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2),
  #                    limits = c(0, 2)) +
  # scale_x_discrete(position = "top") +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "none",
    # legend.title = element_blank(),
    # legend.justification=c(1,0),
    # legend.position=c(0.23, 0.70),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  )
pSNV_prop

# pout <- grid.arrange(p1, p2, nrow = 1, widths = c(6, 6))
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_ann_CpG_SNV_prop.jpg", plot = pout, height = 6, width = 10)

#####

# # calculate 1st 2nd and 3rd quartiles for difference in proportions by annotation
# tmp1 <- aggregate(dt$p_var_dif, list(dt$species, dt$annotation), median)
# colnames(tmp1) <- c("Species", "Category", "Q2")
# tmp2 <- aggregate(dt$p_var_dif, list(dt$species, dt$annotation), quantile, probs = 0.25)
# colnames(tmp2) <- c("Species", "Category", "Q1")
# tmp3 <- aggregate(dt$p_var_dif, list(dt$species, dt$annotation), quantile, probs = 0.75)
# colnames(tmp3) <- c("Species", "Category", "Q3")
# q_dif <- merge(tmp2, tmp1)
# q_dif <- merge(q_dif, tmp3)
# rm(tmp1, tmp2, tmp3)
