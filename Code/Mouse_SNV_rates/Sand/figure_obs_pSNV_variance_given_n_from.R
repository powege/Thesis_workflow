rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)

### IMPORT
dt <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz")

# mean mu rates for k1CG model 
dt <- unique(dt[,c("k1CG_from", "to", "k1CG_mu_rate")])
# simulate
mu_rates <- dt$k1CG_mu_rate
# mu_rates <- 0.001
sim_list <- list()
for(i in 1:length(mu_rates)){
  n_trials <- c(10:10000)
  prob <- rep(mu_rates[i], length(n_trials))
  mean <- n_trials * prob
  variance <- n_trials * prob * ( 1 - prob )
  sd <- sqrt(variance)
  # standardise distance from mean
  sd_standardised <- (sd/mean) 
  sd_min <- (mean/n_trials) - (sd/n_trials)
  sd_max <- (mean/n_trials) + (sd/n_trials)
  sim_list[[i]] <- data.table(prob, n_trials, mean, variance, sd, sd_standardised, sd_min, sd_max)
}
dt_sim <- do.call("rbind", sim_list)
dt_sim$prob <- as.character(dt_sim$prob)

ggplot(dt_sim, aes(x=n_trials, y=sd_standardised, color=prob)) +
  geom_point() +
  ylim(0,2.5) + 
  xlim(0, 500)

min_trials <- rep(NA, length(mu_rates))
min_list <- list()
for(i in 1:length(mu_rates)){
  min_trials[i] <- min(dt_sim$n_trials[dt_sim$prob == mu_rates[i] & dt_sim$sd_min > 0])
  min_list[[i]] <- dt_sim[dt_sim$prob == mu_rates[i] & dt_sim$n_trials == min_trials[i]]
}
dt_min_trials <- do.call("rbind", min_list)




tmp <- dt_sim[n_trials == 100]




