### plot AE sd as opposed to MAE confidennce interval

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(cowplot)

### FUNCTIONS

caterpie <- function(sub){

  k1_from = rep(sub$k1_from[1], 5)
  mutation = rep(sub$mutation[1], 5)
  model = c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer")
  MAE =  c(mean(sub$k1CG_AE),
           mean(sub$k3_AE),
           mean(sub$k5_AE),
           mean(sub$k7_AE),
           mean(sub$k9_AE))
  MAE_CI95 = c(1.96 * (sd(sub$k1CG_AE) / sqrt(length(sub$k1CG_AE))),
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

### SET VARIABLES

snv.file.meth <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_model_fit_specific_CGmethylated.csv"
figure.CGmeth.k1CG.rel <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_k9_CGmeth_k1CGrel_MAE.jpg"

snv.file.unmeth <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_model_fit_specific_CGunmethylated.csv"
figure.CGunmeth.k1CG.rel <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_k9_CGunmeth_k1CGrel_MAE.jpg"


### IMPORT
dt_meth <- fread(snv.file.meth)
dt_unmeth <- fread(snv.file.unmeth)

### FORMAT 

alakazam <- function(dt_meth, status){
dt_meth$mutation <- paste(dt_meth$k1_from, ">", dt_meth$to, sep = " ")
dt_meth_none <- dt_meth[meth_status == "none"]
dt_meth_meth <- dt_meth[meth_status == status]

dt_meth_list <- list()
for (i in 1:length(unique(dt_meth_meth$mutation))){
  sub <- dt_meth_meth[mutation == unique(dt_meth_meth$mutation)[i]]
  dt_meth_list[[i]] <- caterpie(sub)
}
dt_meth_meth <- do.call("rbind", dt_meth_list)
tmp <- ddply(dt_meth_meth, "mutation", metapod_k1CGrel)
dt_meth_meth <- merge(tmp, dt_meth_meth)

dt_meth_none_list <- list()
for (i in 1:length(unique(dt_meth_none$mutation))){
  sub <- dt_meth_none[mutation == unique(dt_meth_none$mutation)[i]]
  dt_meth_none_list[[i]] <- caterpie(sub)
}
dt_meth_none <- do.call("rbind", dt_meth_none_list)
tmp <- ddply(dt_meth_none, "mutation", metapod_k1CGrel)
dt_meth_none <- merge(tmp, dt_meth_none)
rm(tmp)

dt_meth_meth$model <- factor(dt_meth_meth$model, ordered = TRUE, levels = c("1-mer", "A-CG", "3-mer", "5-mer", "7-mer", "9-mer"))
dt_meth_none$model <- factor(dt_meth_none$model, ordered = TRUE, levels = c("1-mer", "A-CG", "3-mer", "5-mer", "7-mer", "9-mer"))
dt_meth_meth <- dt_meth_meth[dt_meth_meth$model %in% c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"),]
dt_meth_none <- dt_meth_none[dt_meth_none$model %in% c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"),]

return(list(dt_meth_meth, dt_meth_none))
}

meth_list <- alakazam(dt_meth, "methylated")
unmeth_list <- alakazam(dt_unmeth, "unmethylated")
dt_meth_meth <- meth_list[[1]]
dt_meth_none <- meth_list[[2]]
dt_unmeth_meth <- unmeth_list[[1]]
dt_unmeth_none <- unmeth_list[[2]]

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
          legend.position = c(0.15, 0.15),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linetype="solid", colour ="black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
}

k1CG_CGmethylated_mae_rel <- plot.mae.k1CGrel(dt_meth_meth, "CG - methylated")
k1CG_CGunmethylated_null_mae_rel <- plot.mae.k1CGrel(dt_meth_none, "CG - null")
out_mae_CGmethylated_rel <- plot_grid(k1CG_CGmethylated_mae_rel, k1CG_CGnull_mae_rel, ncol = 2, nrow = 1)
save_plot(figure.CGmeth.k1CG.rel,
          out_mae_CGmethylated_rel,
          ncol = 1, nrow = 1,
          base_height = 5, base_width = 10)

k1CG_CGunmethylated_mae_rel <- plot.mae.k1CGrel(dt_unmeth_meth, "CG - unmethylated")
k1CG_CGunmethylated_null_mae_rel <- plot.mae.k1CGrel(dt_unmeth_none, "CG - null")
out_mae_CGunmethylated_rel <- plot_grid(k1CG_CGmethylated_mae_rel, k1CG_CGnull_mae_rel, ncol = 2, nrow = 1)
save_plot(figure.CGmeth.k1CG.rel,
          out_mae_CGmethylated_rel,
          ncol = 1, nrow = 1,
          base_height = 5, base_width = 10)


