rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)


### FUNCTIONS
alakazam <- function(dt, lessX){
  
k11 <- dt[,c("k11_from", "to", "k11_from_N", "k11_to_N", "k11_mu_rate")]
k9 <- unique(dt[,c("k9_from", "to", "k9_from_N", "k9_to_N", "k9_mu_rate")])
k7 <- unique(dt[,c("k7_from", "to", "k7_from_N", "k7_to_N", "k7_mu_rate")])
k5 <- unique(dt[,c("k5_from", "to", "k5_from_N", "k5_to_N", "k5_mu_rate")])
k3 <- unique(dt[,c("k3_from", "to", "k3_from_N", "k3_to_N", "k3_mu_rate")])

colnames(k11) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k9) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k7) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k5) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k3) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")

k_list <- list(k3, k5, k7, k9, k11)
total_k <- rep(NA, length(k_list))
n_k_lessX_obs <- rep(NA, length(k_list))
# total_k_changes <- rep(NA, length(k_list))
# n_k_changes_less10_obs <- rep(NA, length(k_list))

for (i in 1:length(k_list)){
  tmp <- k_list[[i]]
  
  total_k[i] <- (nrow(tmp)/3) * 2
  n_k_lessX_obs[i] <- (nrow(tmp[k_from_N < lessX])/3) * 2
  # total_k_changes[i] <- nrow(tmp) * 2
  # n_k_changes_less10_obs[i] <- nrow(tmp[k_to_N < 10]) * 2
  
  print(i)
}

output <- data.table(model = c("3-mer", "5-mer", "7-mer", "9-mer", "11-mer"),
                     total_k,
                     n_k_lessX_obs
                     # total_k_changes,
                     # n_k_changes_less10_obs
                     )

output$prop_k_lessX_obs <- (output$n_k_lessX_obs / output$total_k) * 100
# output$prop_k_changes_less10_obs <- (output$n_k_changes_less10_obs / output$total_k_changes) * 100
output$xaxis_1 <- paste0(output$total_k, "\n(", output$model, ")")
# output$xaxis_2 <- paste0(output$total_k_changes, "\n(", output$model, ")")
output$xaxis_1 <- factor(output$xaxis_1, levels = output$xaxis_1)
# output$xaxis_2 <- factor(output$xaxis_2, levels = output$xaxis_2)

return(output)
}

### IMPORT
dt <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz")
dt_meth <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_methylation_kmer_pSNV_specific.csv.gz")

### FORMAT

lessX <- 100
out_all <- alakazam(dt, lessX)
out_CGmeth <- alakazam(dt_meth[meth_status == "methylated"], lessX)
out_CGunmeth <- alakazam(dt_meth[meth_status == "unmethylated"], lessX)
out_CGmeth$status <- "methylated"
out_CGunmeth$status <- "unmethylated"
out_meth <- rbind(out_CGmeth, out_CGunmeth)

### PLOT 

p_all_k <- ggplot(out_all, aes(x=xaxis_1, y=prop_k_lessX_obs, group = 1)) +
  geom_point() +
  geom_line() +
  xlab("kmers (n)") +
  ylab(paste0("kmers < ", lessX, " observations (%)")) +
  ggtitle("All kmers") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )

# p_all_k_change <- ggplot(out_all, aes(x=xaxis_2, y=prop_k_changes_less10_obs, group = 1)) +
#   geom_point() +
#   geom_line() +
#   xlab("kmer changes (n)") +
#   ylab("kmer changes < 10 observation (%)") +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
#     text = element_text(size=14)
#   )

p_all_k_CG <- ggplot(out_meth, aes(x=xaxis_1, y=prop_k_lessX_obs, group=status)) +
  geom_point(size=2.5) +
  geom_line(aes(linetype=status), size=1.25) +
  xlab("kmers (n)") +
  ylab(paste0("kmers < ", lessX, " observations (%)")) +
  ggtitle("CG kmers") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.2, 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )

# p_all_k_change_CG <- ggplot(out_meth, aes(x=xaxis_2, y=prop_k_changes_less10_obs, group=status)) +
#   geom_point(size=2.5) +
#   geom_line(aes(linetype=status), size=1.25) +
#   xlab("kmer changes (n)") +
#   ylab("kmer changes < 10 observation (%)") +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "top",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
#     text = element_text(size=14)
#   )


pout <- grid.arrange(p_all_k, p_all_k_CG, ncol = 1, nrow = 2)
# pout_all <- grid.arrange(p_all_k, p_all_k_change, ncol = 1, nrow = 2)
# pout_meth <- grid.arrange(p_all_k_CG, p_all_k_change_CG, ncol = 1, nrow = 2)

ggsave("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_kmer_occurence_by_n_parameters.jpg", 
       plot = pout, height = 8, width = 5)

# ggsave("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_kmer_occurence_by_n_parameters.jpg", 
#        plot = pout_all, height = 8, width = 5)
# 
# ggsave("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/figure_kmer_occurence_by_n_parameters_CG_methylation.jpg", 
#        plot = pout_meth, height = 8, width = 5)










