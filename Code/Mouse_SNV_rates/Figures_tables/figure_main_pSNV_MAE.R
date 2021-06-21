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

### SET VARIABLES

AE.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_obs/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_AE.csv"
k3.psnv.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_pSNV.csv.gz"
k1CG.psnv.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_pSNV.csv.gz"
k1.psnv.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1_pSNV.csv.gz"
figure.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_MAE.jpg"

### IMPORT
dt <- fread(AE.infile)
k1_psnv <- fread(paste0("gunzip -cq ", k1.psnv.infile))
k1CG_psnv <- fread(paste0("gunzip -cq ", k1CG.psnv.infile))
k3_psnv <- fread(paste0("gunzip -cq ", k3.psnv.infile))

### FORMAT 
k1CG_psnv <- unique(k1CG_psnv[k1CG_group == "C (CG)" & k1_to == "T"][,c("k1_from", "k1_to", "k1CG_from_N", "k1CG_to_N", "k1CG_mu_rate")])
k1CG_psnv$mutation <- "CG>TG"
colnames(k1CG_psnv) <- c("k1_from", "k1_to", "k1_from_N", "k1_to_N", "k1_mu_rate", "mutation")
k1_psnv$mutation <- paste0(k1_psnv$k1_from, ">", k1_psnv$k1_to)
psnv <- rbind(k1_psnv, k1CG_psnv)
psnv <- psnv[k1_to != "*"]

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

### PLOT BAR
p_psnv <- ggplot(psnv, aes(x=mutation, y=k1_mu_rate, fill=mutation)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = rev(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442"))) +
  xlab("") +
  ylab("Substitution rate") +
  ggtitle("(A)") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p_psnv

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
    ylab("MAE per kb\n(% change relative to 1-mer model)") +
    ggtitle("(B) A>") +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title =  element_text(size = 20, face = "bold"),
          legend.position = c(0.18, 0.25),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linetype="solid", colour ="black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
k1_A_mae_rel

k1_C_mae_rel <- ggplot(dt_C, aes(x=model, y=MAE_k1_rel, group=mutation, color=mutation)) + 
  geom_errorbar(aes(ymax = MAE_k1_rel+MAE_k1_rel_CI95, ymin = MAE_k1_rel-MAE_k1_rel_CI95), 
                stat = "identity",
                position = position_dodge(width = 0.5),
                width=0.2,
                size=1.2) + 
  geom_point(size = 1,
             position = position_dodge(width = 0.5)) + 
  geom_line(position = position_dodge(width = 0.5), size = 1) +
  # scale_colour_manual(values = c("black", "#00AFBB", "#E7B800", "#FC4E07")) +
  xlab("") +
  ylab("MAE per kb\n(% change relative to 1-mer model)") +
  ggtitle("(C) C>") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        legend.position = c(0.8, 0.35),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
k1_C_mae_rel

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
  ggtitle("(D) C>") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title =  element_text(size = 20, face = "bold"),
        legend.position = c(0.18, 0.25),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(linetype="solid", colour ="black"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
k1CG_C_mae_rel

k1CG_CG_mae_rel <- ggplot(dt_CG, aes(x=model, y=MAE_k1CG_rel, group=mutation, color=mutation)) + 
    geom_errorbar(aes(ymax = MAE_k1CG_rel+MAE_k1CG_rel_CI95, ymin = MAE_k1CG_rel-MAE_k1CG_rel_CI95), 
                  stat = "identity",
                  position = position_dodge(width = 0.5),
                  width=0.2,
                  size=1.2) + 
    geom_point(size = 1,
               position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5), size = 1) +
    xlab("") +
    ylab("MAE per kb\n(% change relative to A-CG model)") +
    ggtitle("CG>") +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title =  element_text(size = 20, face = "bold"),
          # legend.position = c(0.2, 0.25),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(linetype="solid", colour ="black"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
k1CG_CG_mae_rel


out_plot <- plot_grid(p_psnv,
                      k1_A_mae_rel, 
                      k1_C_mae_rel, 
                      k1CG_C_mae_rel,
                      ncol = 2, nrow = 2)
save_plot(figure.outfile,
          out_plot,
          ncol = 1, nrow = 1,
          base_height = 8, base_width = 9)


#####


# k3 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_pSNV.csv.gz")
# k3_C <- k3[k1_from =="A" & k1_to == "G"]
# hist(k3_C$k3_mu_rate)
