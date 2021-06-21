rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

dt <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz")
dt <- dt[,c("k7_from", "to", "k7_mu_rate")]
dt <- dt[-which(stri_sub(dt$k7_from, 4, -3) == "CG")] # remove CGs 
dt_A <- dt[which(stri_sub(dt$k7_from, 4, -4) == "A"),]
dt_C <- dt[which(stri_sub(dt$k7_from, 4, -4) == "C"),]
dt_A_wide <- dcast(dt_A, k7_from ~ to, value.var="k7_mu_rate")
dt_C_wide <- dcast(dt_C, k7_from ~ to, value.var="k7_mu_rate")

cor(dt_A_wide[,c("C", "G", "T")], method = "spearman")
cor(dt_C_wide[,c("A", "G", "T")], method = "spearman")

plot(dt_A_wide$C, dt_A_wide$G)
plot(dt_A_wide$C, dt_A_wide$T)
plot(dt_C_wide$A, dt_C_wide$G)
plot(dt_C_wide$A, dt_C_wide$T)



#####



percentile <- ecdf(dt_A$k7_mu_rate)
dt_A$k7_mu_rate_percentile <- ceiling(percentile(dt_A$k7_mu_rate)*100)
dt_A_low <- dt_A[k7_mu_rate_percentile <= 10]
length(unique(dt_A_low$k7_from))
table(dt_C_low$to)

percentile <- ecdf(dt_C$k7_mu_rate)
dt_C$k7_mu_rate_percentile <- ceiling(percentile(dt_C$k7_mu_rate)*100)
dt_C_low <- dt_C[k7_mu_rate_percentile <= 10]
length(unique(dt_C_low$k7_from))
table(dt_A_low$to)




