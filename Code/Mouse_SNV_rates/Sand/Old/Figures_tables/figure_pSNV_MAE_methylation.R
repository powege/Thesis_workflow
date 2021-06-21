rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(cowplot)
library(ggplot2)

### FUNCTIONS
caterpie <- function(sub){
  
  mutation = rep(sub$mutation[1], 2)
  model = c("A-CG", "A-CG_meth")
  MAE =  c(mean(sub$AE_null), mean(sub$AE))
  MAE_CI95 = c(1.96 * (sd(sub$AE_null) / sqrt(length(sub$AE_null))),
               1.96 * (sd(sub$AE) / sqrt(length(sub$AE))))
  output <- data.table(mutation, model, MAE, MAE_CI95)
  
  return(output)
}
metapod_k1CGrel <- function(sub){
  data.table(mutation = sub$mutation,
             model = sub$model,
             MAE_k1CG_rel = ((sub$MAE - sub$MAE[sub$model == "A-CG"]) / sub$MAE[sub$model == "A-CG"])*100,
             MAE_k1CG_rel_CI95 = (sub$MAE_CI95 / sub$MAE[sub$model == "A-CG"])*100
  )
}

### SET VARIABLES

MAE.meth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/Old/MAE_obs/WM_Harr_etal_2016_allSPECIES_CG_resampling_by_chr_model_fit_specific_CG_methylated.csv"
MAE.unmeth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/Old/MAE_obs/WM_Harr_etal_2016_allSPECIES_CG_resampling_by_chr_model_fit_specific_CG_unmethylated.csv"
pSNV.all.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/Old/pSNV/WM_Harr_etal_2016_allSPECIES_CG_pSNV_specific.csv"
pSNV.meth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/Old/pSNV/WM_Harr_etal_2016_allSPECIES_CG_methylated_pSNV_specific.csv"
pSNV.unmeth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/Old/pSNV/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_pSNV_specific.csv"
figure.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/Old/Figures_tables/figure_MAE_methylation.jpg"

### IMPORT
dt_meth <- fread(MAE.meth.infile)
dt_unmeth <- fread(MAE.unmeth.infile)
psnv_all <- fread(pSNV.all.infile)
psnv_meth <- fread(pSNV.meth.infile)
psnv_unmeth <- fread(pSNV.unmeth.infile)

### FORMAT 
dt_meth_list <- list()
for (i in 1:length(unique(dt_meth$mutation))){
  sub <- dt_meth[mutation == unique(dt_meth$mutation)[i]]
  dt_meth_list[[i]] <- caterpie(sub)
}
dt_meth <- do.call("rbind", dt_meth_list)
dt_meth[, 3:4] <- lapply(dt_meth[, 3:4], as.numeric) # convert 
for (i in 1:length(unique(dt_meth$mutation))){
  sub <- dt_meth[mutation == unique(dt_meth$mutation)[i]]
  dt_meth_list[[i]] <- metapod_k1CGrel(sub)
}
dt_meth <- do.call("rbind", dt_meth_list)
dt_meth$model[dt_meth$model == "A-CG_meth"] <- "CG+"

dt_unmeth_list <- list()
for (i in 1:length(unique(dt_unmeth$mutation))){
  sub <- dt_unmeth[mutation == unique(dt_unmeth$mutation)[i]]
  dt_unmeth_list[[i]] <- caterpie(sub)
}
dt_unmeth <- do.call("rbind", dt_unmeth_list)
dt_unmeth[, 3:4] <- lapply(dt_unmeth[, 3:4], as.numeric) # convert to numeric
for (i in 1:length(unique(dt_unmeth$mutation))){
  sub <- dt_unmeth[mutation == unique(dt_unmeth$mutation)[i]]
  dt_unmeth_list[[i]] <- metapod_k1CGrel(sub)
}
dt_unmeth <- do.call("rbind", dt_unmeth_list)
dt_unmeth$model[dt_unmeth$model == "A-CG_meth"] <- "CG-"

dt_mae <- rbind(dt_meth, dt_unmeth)
dt_mae <- dt_mae[model != "A-CG"]
dt_mae$model <- factor(dt_mae$model, levels = c("CG+", "CG-"))


psnv_all <- psnv_all[,c("from", "to", "mu_rate")]
psnv_meth <- psnv_meth[,c("from", "to", "mu_rate")]
psnv_unmeth <- psnv_unmeth[,c("from", "to", "mu_rate")]
psnv_all$status <- "null"
psnv_meth$status <- "CG+"
psnv_unmeth$status <- "CG-"
dt_psnv <- rbind(psnv_all, psnv_meth, psnv_unmeth)
# dt_psnv <- rbind(psnv_meth, psnv_unmeth)
dt_psnv$mutation <- paste0(dt_psnv$from, ">", dt_psnv$to)
dt_psnv$status <- factor(dt_psnv$status, levels = c("null", "CG+", "CG-"))


### PLOT 

fig_psnv <- ggplot(data=dt_psnv, aes(x=status, y=mu_rate, fill=mutation)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  xlab("") +
  ylab("mutation rate") +
  # ggtitle() +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 26, face = "bold"),
        legend.position = c(0.85, 0.84),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
fig_psnv

fig_mae <- ggplot(dt_mae, aes(x=model, y=MAE_k1CG_rel, group=mutation, color=mutation)) + 
    geom_errorbar(aes(ymax = MAE_k1CG_rel+MAE_k1CG_rel_CI95, ymin = MAE_k1CG_rel-MAE_k1CG_rel_CI95), 
                  stat = "identity",
                  position = position_dodge(width = 0.5),
                  width=0.2,
                  size=1.2) + 
    geom_point(size = 1,
               position = position_dodge(width = 0.5)) + 
    xlab("") +
    ylab("MAE per kb\n(% change relative to null model)") +
    # ggtitle() +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title =  element_text(size = 26, face = "bold"),
          legend.position = c(0.83, 0.84),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linetype="solid", colour ="black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
fig_mae
  
out_plot <- plot_grid(fig_psnv, fig_mae,
                      ncol = 2, nrow = 1)

### EXPORT
save_plot(figure.outfile,
          out_plot,
          ncol = 1, nrow = 1,
          base_height = 4, base_width = 8)



