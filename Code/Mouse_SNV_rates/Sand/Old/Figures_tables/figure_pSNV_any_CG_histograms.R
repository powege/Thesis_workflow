# Figure 1 -- Histograms of mutation rates for 3-mer, 5-mer and 7-mer models 
# (ie count ~ mutation rate). Highlight A, C (CG) and C (non-CG) represent 1-mer with lines. 

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(cowplot)

### SET VARS

snv.file.any <- "~/Dropbox/PhD/Data/Mu_rates/MGP_v5_allSTRAIN_kmer_PSNV_any.csv" 

### IMPORT 

kmers <- fread(snv.file.any)

### FORMAT 

kmers <- kmers[k1_from == "A" | k1_from == "C"]

kmers$CG[kmers$CG == "CG"] <- "CG mutation"
kmers$CG[kmers$CG == "nonCG"] <- "nonCG mutation"

kmers$cat <- "C (non CG) >"
kmers$cat[kmers$CG == "CG mutation"] <- "C (CG) >"
kmers$cat[kmers$k1_from == "A"] <- "A >"

k7 <- kmers[,c("k7_from", "k7_mu_rate", "CG", "cat")]
k5 <- unique(kmers[,c("k5_from", "k5_mu_rate", "CG", "cat")])
k3 <- unique(kmers[,c("k3_from", "k3_mu_rate", "CG", "cat")])
k1 <- unique(kmers[,c("k1_from", "k1_mu_rate", "cat")])

k1CG <- kmers$k1CG_mu_rate[kmers$cat == "C (CG) >"][1]
k1nonCG <- kmers$k1CG_mu_rate[kmers$cat == "C (non CG) >"][1]
k1A <- kmers$k1CG_mu_rate[kmers$cat == "A >"][1]
k1C <- kmers$k1_mu_rate[kmers$k1_from == "C"][1]

### PLOT

pk7 <- ggplot(k7, aes(x=k7_mu_rate, fill=cat)) + 
  geom_histogram(binwidth=0.01, colour = "black") + 
  # geom_vline(xintercept = k1CG, linetype = "solid", color = "black", size = 1) +
  # geom_vline(xintercept = k1nonCG, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = k1A, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = k1C, linetype = "dashed", color = "black", size = 1) +
  scale_x_continuous(limits = c(0, 0.3), breaks=c(0, 0.1, 0.2, 0.3)) + 
  # scale_color_manual(values=c("#0072B2", "#D55E00")) +
  xlab("SNV rate") +
  ylab("Count") +
  ggtitle("7-mer (n = 16384)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        # legend.position = "top",
        legend.position = c(0.7, 0.65),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
pk7

pk5 <- ggplot(k5, aes(x=k5_mu_rate, fill=cat)) + 
  geom_histogram(binwidth=0.01, colour = "black") +
  # geom_vline(xintercept = k1CG, linetype = "solid", color = "black", size = 1) +
  # geom_vline(xintercept = k1nonCG, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = k1A, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = k1C, linetype = "dashed", color = "black", size = 1) +
  scale_x_continuous(limits = c(0, 0.3), breaks=c(0, 0.1, 0.2, 0.3)) + 
  # scale_color_manual(values=c("#0072B2", "#D55E00")) +
  xlab("SNV rate") +
  ylab("Count") +
  ggtitle("5-mer (n = 1024)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        # legend.position = "top",
        legend.position = c(0.7, 0.65),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
pk5

pk3 <- ggplot(k3, aes(x=k3_mu_rate, fill=cat)) + 
  geom_histogram(binwidth=0.01, colour = "black") +
  # geom_vline(xintercept = k1CG, linetype = "solid", color = "black", size = 1) +
  # geom_vline(xintercept = k1nonCG, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = k1A, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = k1C, linetype = "dashed", color = "black", size = 1) +
  scale_x_continuous(limits = c(0, 0.3), breaks=c(0, 0.1, 0.2, 0.3)) + 
  # scale_color_manual(values=c("#0072B2", "#D55E00")) +
  xlab("SNV rate") +
  ylab("Count") +
  ggtitle("3-mer (n = 64)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        # legend.position = "top",
        legend.position = c(0.7, 0.65),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
pk3

### EXPORT

output1 <- plot_grid(pk7, pk5, pk3, ncol = 1, nrow = 3)
save_plot("~/Dropbox/PhD/Data/Mu_rates/Figures/Kmer_histograms_PSNV_any.jpg", 
          output1, 
          ncol = 1, nrow = 1, 
          base_height = 8, base_width = 5)
