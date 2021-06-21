rm(list = ls())
graphics.off()

library(stringr)
library(ggplot2)

### Set vars ==========

k7.pSNV.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV.csv.gz"
figure.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/figures_tables/figure_supplementary_nCG_mu_rate.jpg"

### Import data ==========
psnv <- fread(paste0("gunzip -cq ", k7.pSNV.file))

### Format ==========

# identify CG dinucleotides
psnv$CG_status <- "nonCG"
psnv$CG_status[which(stri_sub(psnv$k7_from, 4, -3) == "CG")] <- "CG"

# count number of cytosine or guanine bases in 7-mer
psnv$n_CorG <- str_count(psnv$k7_from, "C") + str_count(psnv$k7_from, "G")

### Model ==========

mod_all <- glm(formula = n_CorG ~ k7_mu_rate, data = psnv, family = "poisson")
summary(mod_all)
coef(summary(mod_all))

non_CG <- psnv[CG_status == "nonCG"]
mod_nonCG <- glm(formula = n_CorG ~ k7_mu_rate, data = non_CG, family = "poisson")
coef(summary(mod_nonCG))

### Plot ==========

psnv$n_CorG <- as.character(psnv$n_CorG)
pA <- ggplot(psnv, aes(x=n_CorG, y=k7_mu_rate)) +
  # geom_boxplot(fill = "grey", outlier.shape = NA) +
  geom_boxplot(fill = "grey") +
  ggtitle("A") +
  xlab("number of C or G in 7-mer") +
  ylab("7-mer substitution rate") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
  )
pA

non_CG$n_CorG <- as.character(non_CG$n_CorG)
pB <- ggplot(non_CG, aes(x=n_CorG, y=k7_mu_rate)) +
  # geom_boxplot(fill = "grey", outlier.shape = NA) +
  geom_boxplot(fill = "grey") +
  ggtitle("B") +
  xlab("number of C or G in 7-mer") +
  ylab("7-mer substitution rate") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
  )
pB

out_plot <- plot_grid(pA,
                      pB,
                      ncol = 1, nrow = 2)

### Export ==========

save_plot(figure.outfile,
          out_plot,
          ncol = 1, nrow = 1,
          base_height = 8, base_width = 5)



##########
# ggplot(psnv, aes(x=n_CG, y=k7_mu_rate)) +
#   geom_boxplot(notch=TRUE)
# sub <- psnv[CG_status != "CG"]
# ggplot(sub, aes(x=n_CG, y=k7_mu_rate)) +
#   geom_boxplot(notch=TRUE)
