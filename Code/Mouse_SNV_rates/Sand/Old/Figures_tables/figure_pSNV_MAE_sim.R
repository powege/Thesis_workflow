### plot AE sd as opposed to MAE confidennce interval

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(cowplot)
library(ggplot2)

### FUNCTIONS

caterpie <- function(sub){
  sub$k1_ae <- abs(sub$exp_mu_k1 - sub$obs_mu)
  sub$k1CG_ae <- abs(sub$exp_mu_k1CG - sub$obs_mu)
  sub$k3_ae <- abs(sub$exp_mu_k3 - sub$obs_mu)
  sub$k5_ae <- abs(sub$exp_mu_k5 - sub$obs_mu)
  sub$k7_ae <- abs(sub$exp_mu_k7 - sub$obs_mu)
  sub$k9_ae <- abs(sub$exp_mu_k9 - sub$obs_mu)

  k1_from = rep(sub$k1_from[1], 6)
  mutation = rep(sub$mutation[1], 6)
  model = c("1-mer", "A-CG", "3-mer", "5-mer", "7-mer", "9-mer")
  MAE =  c(mean(sub$k1_AE),
           mean(sub$k1CG_AE),
           mean(sub$k3_AE),
           mean(sub$k5_AE),
           mean(sub$k7_AE),
           mean(sub$k9_AE))
  MAE_CI95 = c(1.96 * (sd(sub$k1_AE) / sqrt(length(sub$k1_AE))),
               1.96 * (sd(sub$k1CG_AE) / sqrt(length(sub$k1CG_AE))),
               1.96 * (sd(sub$k3_AE) / sqrt(length(sub$k3_AE))),
               1.96 * (sd(sub$k5_AE) / sqrt(length(sub$k5_AE))),
               1.96 * (sd(sub$k7_AE) / sqrt(length(sub$k7_AE))),
               1.96 * (sd(sub$k9_AE) / sqrt(length(sub$k9_AE))))

  output <- data.table( k1_from, mutation, model, MAE, MAE_CI95)
  return(output)
}

metapod_k1CGrel <- function(sub){
  data.table(mutation = sub$mutation,
             model = sub$model,
             MAE_k1CG_rel = ((sub$MAE - sub$MAE[sub$model == "A-CG"]) / sub$MAE[sub$model == "A-CG"])*100,
             MAE_k1CG_rel_CI95 = (sub$MAE_CI95 / sub$MAE[sub$model == "A-CG"])*100
  )
}

beedrill <- function(snv.file.AC, snv.file.CCG, figure.outfile){

### IMPORT
dt_AC <- fread(snv.file.AC)
dt_CCG <- fread(snv.file.CCG)

### FORMAT 
dt_AC$mutation <- paste(dt_AC$k1_from, ">", dt_AC$to, sep = " ")
dt_AC_list <- list()
for (i in 1:length(unique(dt_AC$mutation))){
  sub <- dt_AC[mutation == unique(dt_AC$mutation)[i]]
  dt_AC_list[[i]] <- caterpie(sub)
}
dt_AC <- do.call("rbind", dt_AC_list)
dt_AC[, 4:5] <- lapply(dt_AC[, 4:5], as.numeric) # convert to numeric
tmp <- ddply(dt_AC, "mutation", metapod_k1CGrel)
dt_AC <- merge(tmp, dt_AC)
rm(tmp)
dt_AC <- dt_AC[dt_AC$model %in% c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"),]
dt_AC$model <- factor(dt_AC$model, ordered = TRUE, levels = c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"))
dt_A <- dt_AC[dt_AC$k1_from == "A",]
dt_C <- dt_AC[dt_AC$k1_from == "C",]

dt_CCG$mutation <- paste(dt_CCG$k1_from, ">", dt_CCG$to, sep = " ")
dt_CG <- dt_CCG[CG_status == "CG"]
dt_nonCG <- dt_CCG[CG_status == "nonCG"]
dt_CG_list <- list()
for (i in 1:length(unique(dt_CG$mutation))){
  sub <- dt_CG[mutation == unique(dt_CG$mutation)[i]]
  dt_CG_list[[i]] <- caterpie(sub)
}
dt_CG <- do.call("rbind", dt_CG_list)
tmp <- ddply(dt_CG, "mutation", metapod_k1CGrel)
dt_CG <- merge(tmp, dt_CG)
dt_nonCG_list <- list()
for (i in 1:length(unique(dt_nonCG$mutation))){
  sub <- dt_nonCG[mutation == unique(dt_nonCG$mutation)[i]]
  dt_nonCG_list[[i]] <- caterpie(sub)
}
dt_nonCG <- do.call("rbind", dt_nonCG_list)
tmp <- ddply(dt_nonCG, "mutation", metapod_k1CGrel)
dt_nonCG <- merge(tmp, dt_nonCG)
rm(tmp)
dt_CG$model <- factor(dt_CG$model, ordered = TRUE, levels = c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"))
dt_nonCG$model <- factor(dt_nonCG$model, ordered = TRUE, levels = c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"))
dt_CG <- dt_CG[dt_CG$model %in% c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"),]
dt_nonCG <- dt_nonCG[dt_nonCG$model %in% c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"),]

### PLOT MEAN ABSOLUTE ERROR

plot.mae.k1CGrel <- function(dt_CG, title){
  ggplot(dt_CG, aes(x=model, y=MAE_k1CG_rel, group=mutation, color=mutation)) + 
    geom_errorbar(aes(ymax = MAE_k1CG_rel+MAE_k1CG_rel_CI95, ymin = MAE_k1CG_rel-MAE_k1CG_rel_CI95), 
                  stat = "identity",
                  position = position_dodge(width = 0.5),
                  width=0.2,
                  size=1.2) + 
    geom_point(size = 1,
               position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    xlab("") +
    ylab("MAE per kb\n(% change relative to A-CG model)") +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title =  element_text(size = 26, face = "bold"),
          legend.position = c(0.2, 0.2),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linetype="solid", colour ="black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
}

k1CG_CG_mae_rel <- plot.mae.k1CGrel(dt_CG, "C (CG)")
k1CG_nonCG_mae_rel <- plot.mae.k1CGrel(dt_nonCG, "C (non-CG)")
k1CG_A_mae_rel <- plot.mae.k1CGrel(dt_A, "A")
out_plot <- plot_grid(k1CG_A_mae_rel, k1CG_nonCG_mae_rel, k1CG_CG_mae_rel, ncol = 3, nrow = 1)
save_plot(figure.outfile,
          out_plot,
          ncol = 1, nrow = 1,
          base_height = 5, base_width = 12)

k1CG_A_mae_rel
k1CG_CG_mae_rel
k1CG_nonCG_mae_rel

}

### RUN

# k9 null
beedrill(
snv.file.AC = "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_sim/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_model_fit_specific_AC_k9_null.csv",
snv.file.CCG = "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_sim/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_model_fit_specific_CCG_k9_null.csv",
figure.outfile = "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_MAE_k9_null.jpg"
)

# k3 null
beedrill(
  snv.file.AC = "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_sim/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_model_fit_specific_AC_k3_null.csv",
  snv.file.CCG = "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_sim/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_model_fit_specific_CCG_k3_null.csv",
  figure.outfile = "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_MAE_k3_null.jpg"
)



