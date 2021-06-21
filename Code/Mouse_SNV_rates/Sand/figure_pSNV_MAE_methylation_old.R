rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(cowplot)
library(ggplot2)

caterpie <- function(sub){
  
  k1CG_from = rep(sub$k1CG_from[1], 2)
  mutation = rep(sub$mutation[1], 2)
  model = c("A-CG", "A-CG_meth")
  MAE =  c(mean(sub$k1CG_AE), mean(sub$k1CG_meth_AE))
  MAE_CI95 = c(1.96 * (sd(sub$k1CG_AE) / sqrt(length(sub$k1CG_AE))),
               1.96 * (sd(sub$k1CG_meth_AE) / sqrt(length(sub$k1CG_meth_AE))))
  output <- data.table(k1CG_from, mutation, model, MAE, MAE_CI95)
  
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

MAE.meth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_obs/WM_Harr_etal_2016_allSPECIES_k1CG_resampling_by_chr_model_fit_specific_CG_methylated.csv"
MAE.unmeth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_obs/WM_Harr_etal_2016_allSPECIES_k1CG_resampling_by_chr_model_fit_specific_CG_unmethylated.csv"
pSNV.all.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_pSNV_specific.csv.gz"
pSNV.meth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_CG_methylated_pSNV_specific.csv.gz"
pSNV.unmeth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_CG_unmethylated_pSNV_specific.csv.gz"
figure.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_MAE_methylation.jpg"

### IMPORT
dt_meth <- fread(MAE.meth.infile)
dt_unmeth <- fread(MAE.unmeth.infile)
psnv_all <- fread(paste0("gunzip -cq ",pSNV.all.infile))
psnv_meth <- fread(paste0("gunzip -cq ",pSNV.meth.infile))
psnv_unmeth <- fread(paste0("gunzip -cq ",pSNV.unmeth.infile))

### FORMAT 
dt_meth$mutation <- paste0(dt_meth$k1CG_from, ">", dt_meth$to)
dt_meth_list <- list()
for (i in 1:length(unique(dt_meth$mutation))){
  sub <- dt_meth[mutation == unique(dt_meth$mutation)[i]]
  dt_meth_list[[i]] <- caterpie(sub)
}
dt_meth <- do.call("rbind", dt_meth_list)
dt_meth[, 4:5] <- lapply(dt_meth[, 4:5], as.numeric) # convert to numeric
dt_meth_list <- list()
for (i in 1:length(unique(dt_meth$mutation))){
  sub <- dt_meth[mutation == unique(dt_meth$mutation)[i]]
  dt_meth_list[[i]] <- metapod_k1CGrel(sub)
}
dt_meth <- do.call("rbind", dt_meth_list)
dt_meth$model[dt_meth$model == "A-CG_meth"] <- "CG+"

dt_unmeth$mutation <- paste0(dt_unmeth$k1CG_from, ">", dt_unmeth$to)
dt_unmeth_list <- list()
for (i in 1:length(unique(dt_unmeth$mutation))){
  sub <- dt_unmeth[mutation == unique(dt_unmeth$mutation)[i]]
  dt_unmeth_list[[i]] <- caterpie(sub)
}
dt_unmeth <- do.call("rbind", dt_unmeth_list)
dt_unmeth[, 4:5] <- lapply(dt_unmeth[, 4:5], as.numeric) # convert to numeric
dt_unmeth_list <- list()
for (i in 1:length(unique(dt_unmeth$mutation))){
  sub <- dt_unmeth[mutation == unique(dt_unmeth$mutation)[i]]
  dt_unmeth_list[[i]] <- metapod_k1CGrel(sub)
}
dt_unmeth <- do.call("rbind", dt_unmeth_list)
dt_unmeth$model[dt_unmeth$model == "A-CG_meth"] <- "CG-"

dt_mae <- rbind(dt_meth, dt_unmeth)
dt_mae <- dt_mae[model != "A-CG"]

psnv_all <- psnv_all[k1CG_from == "C (CG)"]
psnv_all$k1CG_from <- "C"
psnv_all <- unique(psnv_all[,c("k1CG_from", "to", "k1CG_mu_rate")])
psnv_meth <- psnv_meth[,c("k1CG_from", "to", "k1CG_mu_rate")]
psnv_unmeth <- psnv_unmeth[,c("k1CG_from", "to", "k1CG_mu_rate")]
psnv_all$status <- "all"
psnv_meth$status <- "CG+"
psnv_unmeth$status <- "CG-"
dt_psnv <- rbind(psnv_all, psnv_meth, psnv_unmeth)
# dt_psnv <- rbind(psnv_meth, psnv_unmeth)
dt_psnv$mutation <- paste0(dt_psnv$k1CG_from, ">", dt_psnv$to)


fig_psnv <- ggplot(data=dt_psnv, aes(x=status, y=k1CG_mu_rate, fill=mutation)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  xlab("") +
  ylab("mutation rate") +
  # ggtitle() +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 26, face = "bold"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


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
          legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linetype="solid", colour ="black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )

  
out_plot <- plot_grid(fig_psnv, fig_mae,
                      ncol = 2, nrow = 1)
save_plot(figure.outfile,
          out_plot,
          ncol = 1, nrow = 1,
          base_height = 4, base_width = 8)



#####

#####

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
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
fig_mae







s
  
