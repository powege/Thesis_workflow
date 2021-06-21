rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS
alakazam <- function(dt){
  
k11 <- dt[,c("k11_from", "to", "k11_from_N", "k11_to_N", "k11_mu_rate")]
k9 <- unique(dt[,c("k9_from", "to", "k9_from_N", "k9_to_N", "k9_mu_rate")])
k7 <- unique(dt[,c("k7_from", "to", "k7_from_N", "k7_to_N", "k7_mu_rate")])
k5 <- unique(dt[,c("k5_from", "to", "k5_from_N", "k5_to_N", "k5_mu_rate")])
k3 <- unique(dt[,c("k3_from", "to", "k3_from_N", "k3_to_N", "k3_mu_rate")])
k1CG <- unique(dt[,c("k1CG_from", "to", "k1CG_from_N", "k1CG_to_N", "k1CG_mu_rate")])

colnames(k11) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k9) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k7) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k5) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k3) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
colnames(k1CG) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")


k_list <- list(k11, k9, k7, k5, k3, k1CG)

total_k <- rep(NA, length(k_list))
total_k_changes <- rep(NA, length(k_list))

k_from_mean <- rep(NA, length(k_list))
k_from_median <- rep(NA, length(k_list))
k_from_min <- rep(NA, length(k_list))
k_from_max <- rep(NA, length(k_list))
k_from_total <- rep(NA, length(k_list))

k_to_mean <- rep(NA, length(k_list))
k_to_median <- rep(NA, length(k_list))
k_to_min <- rep(NA, length(k_list))
k_to_max <- rep(NA, length(k_list))
k_to_total <- rep(NA, length(k_list))

mu_rate_mean <- rep(NA, length(k_list))
mu_rate_median <- rep(NA, length(k_list))
mu_rate_min <- rep(NA, length(k_list))
mu_rate_max <- rep(NA, length(k_list))

for (i in 1:length(k_list)){
  tmp <- k_list[[i]]
  
  total_k[i] <- (nrow(tmp) / 3) * 2
  total_k_changes[i] <- nrow(tmp) * 2
  
  k_from_mean[i] <- mean(tmp$k_from_N)
  k_from_median[i] <- median(tmp$k_from_N)
  k_from_min[i] <- min(tmp$k_from_N)
  k_from_max[i] <- max(tmp$k_from_N)
  k_from_total[i] <- sum(tmp$k_from_N) / 3
  
  k_to_mean[i] <- mean(tmp$k_to_N)
  k_to_median[i] <- median(tmp$k_to_N)
  k_to_min[i] <- min(tmp$k_to_N)
  k_to_max[i] <- max(tmp$k_to_N)
  k_to_total[i] <- sum(tmp$k_to_N)
  
  mu_rate_mean[i] <- mean(tmp$k_mu_rate)
  mu_rate_median[i] <- median(tmp$k_mu_rate)
  mu_rate_min[i] <- min(tmp$k_mu_rate)
  mu_rate_max[i] <- max(tmp$k_mu_rate)
  
  print(i)
}

output <- data.table(model = c("k11", "k9", "k7", "k5", "k3", "k1CG"),
                     total_k,
                     total_k_changes,
                     k_from_mean,
                     k_from_median,
                     k_from_min,
                     k_from_max,
                     k_from_total,
                     k_to_mean,
                     k_to_median,
                     k_to_min,
                     k_to_max,
                     k_to_total,
                     mu_rate_mean,
                     mu_rate_median,
                     mu_rate_min,
                     mu_rate_max)

return(output)
}

### IMPORT
dt <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz")
dt_meth <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_methylation_kmer_pSNV_specific.csv.gz")

### FORMAT
out_all <- alakazam(dt)
out_CGmeth <- alakazam(dt_meth[meth_status == "methylated"])
out_CGunmeth <- alakazam(dt_meth[meth_status == "unmethylated"])

out_all$meth_status <- "none"
out_CGmeth$meth_status <- "methylated"
out_CGunmeth$meth_status <- "unmethylated"

out_table <- rbind(out_all, out_CGmeth, out_CGunmeth)

### EXPORT
fwrite(out_table, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/table_kmer_SNV_counts_summary.csv")

















