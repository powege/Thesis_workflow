# SCRIPT that plots:
#   the proportion of CG dinucleotides by annotation
#   the proportion of SNVs by annotation 

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)

### SET VARS
order <-  rev(c("Exon - CDS", 
                "Exon - UTR",
                "Exon - other",
                "Promoter", 
                "Enhancer - proximal", 
                "Enhancer - distal", 
                "CTCF binding", 
                "Miscellaneous", 
                "Intron - proximal",
                "Intron - distal",
                "Unannotated"))

### IMPORT
  mvar <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_GRC38_GENCODE_RegBuild_annotation_SNV_CG_total.csv")
  hvar <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Human_GRC38_GENCODE_RegBuild_annotation_SNV_CG_total.csv")
  
  ### FORMAT
  
  # colnames
  colnames(hvar) <- c("chromosome", "annotation", "n_bp", "n_CG", "n_SNV")
  colnames(mvar) <- c("chromosome", "annotation", "n_bp", "n_CG", "n_SNV")
  
  # totals by annotation
  # human
  annotation <- unique(hvar$annotation)
  n_bp <- rep(NA, length(annotation))
  n_CG <- rep(NA, length(annotation))
  n_SNV <- rep(NA, length(annotation))
  for (i in 1:length(annotation)){
    tmp_sub <- hvar[annotation == annotation[i]]
    n_bp[i] <- sum(tmp_sub$n_bp)
    n_CG[i] <- sum(tmp_sub$n_CG)
    n_SNV[i] <- sum(tmp_sub$n_SNV)
  }
  hvar <- data.table(annotation, 
                     n_bp,
                     n_CG, 
                     n_SNV)
  # mouse
  annotation <- unique(mvar$annotation)
  n_bp <- rep(NA, length(annotation))
  n_CG <- rep(NA, length(annotation))
  n_SNV <- rep(NA, length(annotation))
  for (i in 1:length(annotation)){
    tmp_sub <- mvar[annotation == annotation[i]]
    n_bp[i] <- sum(tmp_sub$n_bp)
    n_CG[i] <- sum(tmp_sub$n_CG)
    n_SNV[i] <- sum(tmp_sub$n_SNV)
  }
  mvar <- data.table(annotation, 
                     n_bp,
                     n_CG, 
                     n_SNV)
  rm(tmp_sub)
  
  # calculate totals
  h_CG_tot <- (hvar$n_CG[hvar$annotation == "total"] / hvar$n_bp[hvar$annotation == "total"]) * 1000
  h_SNV_tot <-(hvar$n_SNV[hvar$annotation == "total"] / hvar$n_bp[hvar$annotation == "total"]) * 1000
  m_CG_tot <- (mvar$n_CG[mvar$annotation == "total"] / mvar$n_bp[mvar$annotation == "total"]) * 1000
  m_SNV_tot <-(mvar$n_SNV[mvar$annotation == "total"] / mvar$n_bp[mvar$annotation == "total"]) * 1000
  
  # subset annotations
  hvar <-  hvar[annotation != "total"]
  mvar <- mvar[annotation != "total"]
  
  # calculate proportions by sequence
  hvar$p_CG <- (hvar$n_CG / hvar$n_bp) * 1000
  mvar$p_CG <- (mvar$n_CG / mvar$n_bp) * 1000
  hvar$p_SNV <- (hvar$n_SNV / hvar$n_bp) * 1000
  mvar$p_SNV <- (mvar$n_SNV / mvar$n_bp) * 1000
  
  # calculate fold change
  hvar$fc_CG <- hvar$p_CG/h_CG_tot
  hvar$fc_SNV <- hvar$p_SNV/h_SNV_tot
  mvar$fc_CG <- mvar$p_CG/m_CG_tot
  mvar$fc_SNV <- mvar$p_SNV/m_SNV_tot
  
  # add species
  hvar$species <- "Human"
  mvar$species <- "Mouse"
  
  # rbind
  dt <- rbind(hvar, mvar)
  rm(hvar, mvar)
  
  # subset annotations
  dt <- dt[annotation %in% order]
  
  # # set rectangle coordinates 
  # rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1), 
  #                     xmax = tail(seq, -1), 
  #                     annotation = levels(as.factor(order)),
  #                     rect_type = c("a", "c"))
  # rects <- rbind(rects)
  # dt <- merge(dt, rects)

  # set factor order
  dt$annotation <- factor(dt$annotation, levels = order)
  
  # subset mouse
  m_dt <- dt[species == "Mouse"]
  
  
  ### PLOT MOUSE
  
  # CG & SNV
  m_dt_long <- melt(m_dt, id.vars = "annotation", measure.vars = c("fc_CG", "fc_SNV"))
  m_dt_long$variable <- as.character(  m_dt_long$variable)
  m_dt_long$variable[m_dt_long$variable == "fc_CG"] <- "CG dinucleotides\nper kb"
  m_dt_long$variable[m_dt_long$variable == "fc_SNV"] <- "SNVs per kb"
  
  p_CG_SNV_m <- ggplot(m_dt_long, aes(x=annotation, y=value, fill=variable)) + 
    geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
    geom_hline(yintercept=1,
               linetype="dashed",
               color = "black",
               size=1) +
    xlab("Genomic annotation") +
    ylab("Fold change relative to genomic average") +
    # ggtitle("CG dinucleotides") +
    # scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2),
    #                    limits = c(0, 2)) +
    # scale_x_discrete(position = "top") +
    coord_flip() +
    theme_bw() +
    theme(
      # legend.position = "none",
      legend.title = element_blank(),
      # legend.justification=c(1,0),
      # legend.position=c(0.23, 0.70),
      # legend.box.background = element_rect(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size=14)
    )
  p_CG_SNV_m

  # CG
  pCG_m <- ggplot(m_dt, aes(x=annotation, y=p_CG, fill=annotation)) + 
    geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
    geom_hline(yintercept=m_CG_tot,
               linetype="dashed",
               color = "black",
               size=1) +
    xlab("Genomic annotation") +
    ylab("CG dinucleotides per kb") +
    scale_fill_brewer(palette="Set3", direction = -1) +
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
  pCG_m
  
  # SNV
  pSNV_m <- ggplot(m_dt, aes(x=annotation, y=p_SNV, fill=annotation)) + 
    geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
    geom_hline(yintercept=m_SNV_tot,
               linetype="dashed",
               color = "black",
               size=1) +
    xlab("Genomic annotation") +
    ylab("SNV sites per kb") +
    scale_fill_brewer(palette="Set3", direction=-1) +
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
  pSNV_m
  
  pout <- grid.arrange(pSNV_m, pCG_m, nrow = 1, widths = c(5, 5))
  ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_mouse_CG_SNV_totals.jpg", plot = pout, height = 5, width = 10)
  
  
### PLOT COMBINED

# CG
pCG <- ggplot(dt, aes(x=annotation, y=p_CG, fill=species)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
  # geom_rect(
  #   aes(xmin = dt$xmin,
  #       xmax = dt$xmax,
  #       ymin = -Inf,
  #       ymax = Inf,
  #       fill = dt$rect_type),
  #   color = NA,
  #   alpha = 0.5,
  #   show.legend = F) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) + 
  geom_hline(yintercept=h_CG_tot,
             linetype="dashed",
             color = "red",
             size=1) +
  geom_hline(yintercept=m_CG_tot,
             linetype="dashed",
             color = "blue",
             size=1) +
  # scale_fill_manual(values = c("grey", "white", "red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
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
pCG

# pout <- grid.arrange(p1, p2, nrow = 1, widths = c(6, 6))
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_ann_CpG_SNV_prop.jpg", plot = pout, height = 6, width = 10)

# SNV
pSNV <- ggplot(dt, aes(x=annotation, y=p_SNV, fill=species)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  # geom_rect(
  #   aes(xmin = dt$xmin,
  #       xmax = dt$xmax,
  #       ymin = -Inf,
  #       ymax = Inf,
  #       fill = dt$rect_type),
  #   color = NA,
  #   alpha = 0.5,
  #   show.legend = F) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) + 
  geom_hline(yintercept=h_SNV_tot,
             linetype="dashed",
             color = "red",
             size=1) +
  geom_hline(yintercept=m_SNV_tot,
             linetype="dashed",
             color = "blue",
             size=1) +
  # scale_fill_manual(values = c("grey", "white", "red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  xlab("Genomic annotation") +
  ylab("SNV sites per kb") +
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
pSNV

# pout <- grid.arrange(p1, p2, nrow = 1, widths = c(6, 6))
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_ann_CpG_SNV_prop.jpg", plot = pout, height = 6, width = 10)

#####

