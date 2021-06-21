rm(list = ls())
graphics.off()

library(data.table)

### Import data ==========

k9 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k9_pSNV.csv.gz")
k7 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV.csv.gz")
k5 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k5_pSNV.csv.gz")
k3 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_pSNV.csv.gz")
acg <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_pSNV.csv.gz")
k1 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1_pSNV.csv.gz")

### Format ==========

# k3
k3_summary <- k3[, N_to_Q1:=quantile(k3_to_N, 0.25), by=.(k1_from, k1_to)] 
k3_summary <- k3_summary[, N_to_median:=median(k3_to_N), by=.(k1_from, k1_to)] 
k3_summary <- k3_summary[, N_to_Q3:=quantile(k3_to_N, 0.75), by=.(k1_from, k1_to)] 
k3_summary <- k3_summary[, N_to_mean:=mean(k3_to_N), by=.(k1_from, k1_to)] 

k3_summary <- k3[, N_from_Q1:=quantile(k3_from_N, 0.25), by=.(k1_from, k1_to)] 
k3_summary <- k3_summary[, N_from_median:=median(k3_from_N), by=.(k1_from, k1_to)] 
k3_summary <- k3_summary[, N_from_Q3:=quantile(k3_from_N, 0.75), by=.(k1_from, k1_to)] 
k3_summary <- k3_summary[, N_from_mean:=mean(k3_from_N), by=.(k1_from, k1_to)] 

k3_summary <- unique(k3_summary[,c("k1_from","k1_to","N_to_Q1","N_to_median","N_to_Q3","N_to_mean","N_from_Q1","N_from_median","N_from_Q3","N_from_mean")])
k3_summary$kmer_model <- "3-mer"

# k5
k5_summary <- k5[, N_to_Q1:=quantile(k5_to_N, 0.25), by=.(k1_from, k1_to)] 
k5_summary <- k5_summary[, N_to_median:=median(k5_to_N), by=.(k1_from, k1_to)] 
k5_summary <- k5_summary[, N_to_Q3:=quantile(k5_to_N, 0.75), by=.(k1_from, k1_to)] 
k5_summary <- k5_summary[, N_to_mean:=mean(k5_to_N), by=.(k1_from, k1_to)] 

k5_summary <- k5[, N_from_Q1:=quantile(k5_from_N, 0.25), by=.(k1_from, k1_to)] 
k5_summary <- k5_summary[, N_from_median:=median(k5_from_N), by=.(k1_from, k1_to)] 
k5_summary <- k5_summary[, N_from_Q3:=quantile(k5_from_N, 0.75), by=.(k1_from, k1_to)] 
k5_summary <- k5_summary[, N_from_mean:=mean(k5_from_N), by=.(k1_from, k1_to)] 

k5_summary <- unique(k5_summary[,c("k1_from","k1_to","N_to_Q1","N_to_median","N_to_Q3","N_to_mean","N_from_Q1","N_from_median","N_from_Q3","N_from_mean")])
k5_summary$kmer_model <- "5-mer"

# k7
k7_summary <- k7[, N_to_Q1:=quantile(k7_to_N, 0.25), by=.(k1_from, k1_to)] 
k7_summary <- k7_summary[, N_to_median:=median(k7_to_N), by=.(k1_from, k1_to)] 
k7_summary <- k7_summary[, N_to_Q3:=quantile(k7_to_N, 0.75), by=.(k1_from, k1_to)] 
k7_summary <- k7_summary[, N_to_mean:=mean(k7_to_N), by=.(k1_from, k1_to)] 

k7_summary <- k7[, N_from_Q1:=quantile(k7_from_N, 0.25), by=.(k1_from, k1_to)] 
k7_summary <- k7_summary[, N_from_median:=median(k7_from_N), by=.(k1_from, k1_to)] 
k7_summary <- k7_summary[, N_from_Q3:=quantile(k7_from_N, 0.75), by=.(k1_from, k1_to)] 
k7_summary <- k7_summary[, N_from_mean:=mean(k7_from_N), by=.(k1_from, k1_to)] 

k7_summary <- unique(k7_summary[,c("k1_from","k1_to","N_to_Q1","N_to_median","N_to_Q3","N_to_mean","N_from_Q1","N_from_median","N_from_Q3","N_from_mean")])
k7_summary$kmer_model <- "7-mer"

# k9
k9_summary <- k9[, N_to_Q1:=quantile(k9_to_N, 0.25), by=.(k1_from, k1_to)] 
k9_summary <- k9_summary[, N_to_median:=median(k9_to_N), by=.(k1_from, k1_to)] 
k9_summary <- k9_summary[, N_to_Q3:=quantile(k9_to_N, 0.75), by=.(k1_from, k1_to)] 
k9_summary <- k9_summary[, N_to_mean:=mean(k9_to_N), by=.(k1_from, k1_to)] 

k9_summary <- k9[, N_from_Q1:=quantile(k9_from_N, 0.25), by=.(k1_from, k1_to)] 
k9_summary <- k9_summary[, N_from_median:=median(k9_from_N), by=.(k1_from, k1_to)] 
k9_summary <- k9_summary[, N_from_Q3:=quantile(k9_from_N, 0.75), by=.(k1_from, k1_to)] 
k9_summary <- k9_summary[, N_from_mean:=mean(k9_from_N), by=.(k1_from, k1_to)] 

k9_summary <- unique(k9_summary[,c("k1_from","k1_to","N_to_Q1","N_to_median","N_to_Q3","N_to_mean","N_from_Q1","N_from_median","N_from_Q3","N_from_mean")])
k9_summary$kmer_model <- "9-mer"

table_out <- rbind(k3_summary, k5_summary, k7_summary, k9_summary)


### Export data ==========

fwrite(table_out, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/figures_tables/table_supplementary_kmer_count_IQR.csv")





