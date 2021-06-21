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

kakuna_k1rel <- function(sub){
  data.table(mutation = sub$mutation,
             model = sub$model,
             MAE_k1_rel = ((sub$MAE - sub$MAE[sub$model == "1-mer"]) / sub$MAE[sub$model == "1-mer"])*100,
             MAE_k1_rel_CI95 = (sub$MAE_CI95 / sub$MAE[sub$model == "1-mer"])*100
  )
}

metapod_k1CGrel <- function(sub){
  data.table(mutation = sub$mutation,
             model = sub$model,
             MAE_k1CG_rel = ((sub$MAE - sub$MAE[sub$model == "A-CG"]) / sub$MAE[sub$model == "A-CG"])*100,
             MAE_k1CG_rel_CI95 = (sub$MAE_CI95 / sub$MAE[sub$model == "A-CG"])*100
  )
}

alltogether <- function(AE.infile, titleA, titleC){
  
  ### IMPORT
  dt <- fread(AE.infile)
  
  ### FORMAT 
  dt$mutation <- paste(dt$k1_from, ">", dt$k1_to, sep = " ")
  dt_AC <- dt[k1_from %in% c("A", "C")]
  dt_CCG <- dt[k1_from %in% c("C", "CG")]
  
  dt_AC_list <- list()
  for (i in 1:length(unique(dt_AC$mutation))){
    sub <- dt_AC[mutation == unique(dt_AC$mutation)[i]]
    dt_AC_list[[i]] <- caterpie(sub)
  }
  dt_AC <- do.call("rbind", dt_AC_list)
  dt_AC[, 4:5] <- lapply(dt_AC[, 4:5], as.numeric) # convert to numeric
  tmp <- ddply(dt_AC, "mutation", kakuna_k1rel)
  dt_AC <- merge(tmp, dt_AC)
  dt_AC$model <- factor(dt_AC$model, ordered = TRUE, levels = c("1-mer", "A-CG", "3-mer", "5-mer", "7-mer", "9-mer"))
  dt_A <- dt_AC[dt_AC$k1_from == "A",]
  dt_C <- dt_AC[dt_AC$k1_from == "C",]
  
  dt_CCG_list <- list()
  for (i in 1:length(unique(dt_CCG$mutation))){
    sub <- dt_CCG[mutation == unique(dt_CCG$mutation)[i]]
    dt_CCG_list[[i]] <- caterpie(sub)
  }
  dt_CCG <- do.call("rbind", dt_CCG_list)
  tmp <- ddply(dt_CCG, "mutation", metapod_k1CGrel)
  dt_CCG <- merge(tmp, dt_CCG)
  dt_CCG$model <- factor(dt_CCG$model, ordered = TRUE, levels = c("1-mer", "A-CG", "3-mer", "5-mer", "7-mer", "9-mer"))
  dt_CCG <- dt_CCG[dt_CCG$model %in% c("A-CG", "3-mer", "5-mer", "7-mer", "9-mer"),]
  dt_C2 <- dt_CCG[dt_CCG$k1_from == "C",]
  dt_CG <- dt_CCG[dt_CCG$k1_from == "CG",]
  
  dt_A <- dt_A[dt_A$model != "1-mer" & dt_A$mutation != "A > *",]
  dt_C2 <- dt_C2[dt_C2$model != "1-mer" & dt_C2$mutation != "C > *",]
  
  
  ### PLOT MEAN ABSOLUTE ERROR
  k1_A_mae_rel <- ggplot(dt_A, aes(x=model, y=MAE_k1_rel, group=mutation, color=mutation)) + 
    geom_errorbar(aes(ymax = MAE_k1_rel+MAE_k1_rel_CI95, ymin = MAE_k1_rel-MAE_k1_rel_CI95), 
                  stat = "identity",
                  position = position_dodge(width = 0.5),
                  width=0.2,
                  size=1.2) + 
    geom_point(size = 1,
               position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5), size = 1) +
    # scale_colour_manual(values = c("black", "#E69F00", "#56B4E9", "#009E73")) +
    xlab("") +
    ylab("MAE per kb\n(% change relative to A-CG model)") +
    ggtitle(titleA) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title =  element_text(size = 20, face = "bold"),
          legend.position = c(0.15, 0.17),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linetype="solid", colour ="black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  k1_A_mae_rel
  
  # k1_C_mae_rel <- ggplot(dt_C, aes(x=model, y=MAE_k1_rel, group=mutation, color=mutation)) + 
  #   geom_errorbar(aes(ymax = MAE_k1_rel+MAE_k1_rel_CI95, ymin = MAE_k1_rel-MAE_k1_rel_CI95), 
  #                 stat = "identity",
  #                 position = position_dodge(width = 0.5),
  #                 width=0.2,
  #                 size=1.2) + 
  #   geom_point(size = 1,
  #              position = position_dodge(width = 0.5)) + 
  #   geom_line(position = position_dodge(width = 0.5), size = 1) +
  #   # scale_colour_manual(values = c("black", "#00AFBB", "#E7B800", "#FC4E07")) +
  #   xlab("") +
  #   ylab("MAE per kb\n(% change relative to 1-mer model)") +
  #   ggtitle("(B) C>") +
  #   theme_classic() +
  #   theme(axis.text = element_text(size = 14),
  #         axis.title = element_text(size = 14),
  #         plot.title =  element_text(size = 20, face = "bold"),
  #         legend.position = c(0.8, 0.35),
  #         legend.title = element_blank(),
  #         legend.text = element_text(size = 14),
  #         legend.background = element_rect(linetype="solid", colour ="black"),
  #         plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1)
  #   )
  # k1_C_mae_rel
  
  k1CG_C_mae_rel <- ggplot(dt_C2, aes(x=model, y=MAE_k1CG_rel, group=mutation, color=mutation)) + 
    geom_errorbar(aes(ymax = MAE_k1CG_rel+MAE_k1CG_rel_CI95, ymin = MAE_k1CG_rel-MAE_k1CG_rel_CI95), 
                  stat = "identity",
                  position = position_dodge(width = 0.5),
                  width=0.2,
                  size=1.2) + 
    geom_point(size = 1,
               position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5), size = 1) +
    # scale_colour_manual(values = c("black", "#0072B2", "#D55E00", "#CC79A7")) +
    xlab("") +
    ylab("MAE per kb\n(% change relative to A-CG model)") +
    ggtitle(titleC) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title =  element_text(size = 20, face = "bold"),
          legend.position = c(0.15, 0.17),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linetype="solid", colour ="black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  k1CG_C_mae_rel
  
  # k1CG_CG_mae_rel <- ggplot(dt_CG, aes(x=model, y=MAE_k1CG_rel, group=mutation, color=mutation)) + 
  #     geom_errorbar(aes(ymax = MAE_k1CG_rel+MAE_k1CG_rel_CI95, ymin = MAE_k1CG_rel-MAE_k1CG_rel_CI95), 
  #                   stat = "identity",
  #                   position = position_dodge(width = 0.5),
  #                   width=0.2,
  #                   size=1.2) + 
  #     geom_point(size = 1,
  #                position = position_dodge(width = 0.5)) + 
  #     geom_line(position = position_dodge(width = 0.5), size = 1) +
  #     xlab("") +
  #     ylab("MAE per kb\n(% change relative to A-CG model)") +
  #     ggtitle("CG>") +
  #     theme_classic() +
  #     theme(axis.text = element_text(size = 14),
  #           axis.title = element_text(size = 14),
  #           plot.title =  element_text(size = 20, face = "bold"),
  #           # legend.position = c(0.2, 0.25),
  #           legend.title = element_blank(),
  #           legend.text = element_text(size = 14),
  #           legend.background = element_rect(linetype="solid", colour ="black"),
  #           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
  #           panel.border = element_rect(colour = "black", fill=NA, size=1)
  #     )
  # k1CG_CG_mae_rel
  
  return(list(k1_A_mae_rel, k1CG_C_mae_rel))
}

### SET VARIABLES

AE.infile.k3 <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_sim/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_AE_k3_sim.csv"
AE.infile.k9 <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_sim/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_AE_k9_sim.csv"
figure.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_supplement_MAE_sim.jpg"

### RUN

k3_figs <- alltogether(AE.infile.k3, titleA = "(A) 3-mer null", titleC = "(B) 3-mer null")
pan_A <- k3_figs[[1]]
pan_B <- k3_figs[[2]]

k9_figs <- alltogether(AE.infile.k9, titleA = "(C) 9-mer null", titleC = "(D) 9-mer null")
pan_C <- k9_figs[[1]]
pan_D <- k9_figs[[2]]


out_plot <- plot_grid(pan_A,
                      pan_B,
                      pan_C,
                      pan_D,
                      ncol = 2, nrow = 2)
save_plot(figure.outfile,
          out_plot,
          ncol = 1, nrow = 1,
          base_height = 8, base_width = 9)


#####


# k3 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_pSNV.csv.gz")
# k3_C <- k3[k1_from =="A" & k1_to == "G"]
# hist(k3_C$k3_mu_rate)
