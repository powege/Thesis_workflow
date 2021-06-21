### SCRIPT that creates simulated datasets with same k9_from_N as observed:
# 1. where k9_to_N is sampled from a binomial distribution with same pMU as k3 model
# 2. where k9_to_N is sampled from a binomial distribution with same pMU as k9 model 

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 

### SET VARS
k9.counts.chr.infile <- args[1]
k9.p.snv.infile <- args[2]
k3.p.snv.infile <- args[3]
k9.null.outfile <- args[4]
k3.null.outfile <- args[5]

# k9.counts.chr.infile <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k9_SNV_counts_by_chr.csv.gz"
# k3.p.snv.infile <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_pSNV_specific.csv.gz"
# k9.p.snv.infile <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k9_pSNV_specific.csv.gz"
# k3.null.outfile <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k9_SNV_counts_by_chr_k3_sim.csv.gz"
# k9.null.outfile <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k9_SNV_counts_by_chr_k9_sim.csv.gz"

### IMPORT 
k9_chr <- fread(paste0("gunzip -cq ", k9.counts.chr.infile))
k9_p_snv <- fread(paste0("gunzip -cq ", k9.p.snv.infile))
k3_p_snv <- fread(paste0("gunzip -cq ", k3.p.snv.infile))


### FORMAT
k9_chr$k1_to <- stri_sub(k9_chr$k9_to, 5, -5)

k9_p_snv <- k9_p_snv[k1_to %in% c("A", "C", "G", "T")]
k3_p_snv <- k3_p_snv[k1_to %in% c("A", "C", "G", "T")]

k9_p_snv$k3_from <- stri_sub(k9_p_snv$k9_from, 4, -4)
k9_p_snv <- k3_p_snv[k9_p_snv, on = c("k3_from", "k1_to")]

k9 <- unique(k9_p_snv[,c("k9_from", "k1_to", "k9_mu_rate")])
k3 <- unique(k9_p_snv[,c("k9_from", "k1_to", "k3_mu_rate")])

# k9$k9_mu_rate[k9$k9_mu_rate == 0] <- min(k9$k9_mu_rate[k9$k9_mu_rate > 0])
# k3$k3_mu_rate[k3$k3_mu_rate == 0] <- min(k3$k3_mu_rate[k3$k3_mu_rate > 0])

k9_null <- k9_chr[k9, on = c("k9_from", "k1_to")]
k3_null <- k9_chr[k3, on = c("k9_from", "k1_to")]

k9_null$k9_to_N <- rbinom(size = k9_null$k9_from_N, n = nrow(k9_null), prob = k9_null$k9_mu_rate)
k3_null$k9_to_N <- rbinom(size = k3_null$k9_from_N, n = nrow(k3_null), prob = k3_null$k3_mu_rate)

k9_null <- k9_null[,c("chromosome", "k9_from", "k9_to", "k9_from_N", "k9_to_N")]
k3_null <- k3_null[,c("chromosome", "k9_from", "k9_to", "k9_from_N", "k9_to_N")]

### EXPORT 
fwrite(k9_null, k9.null.outfile, compress = "gzip")
fwrite(k3_null, k3.null.outfile, compress = "gzip")






