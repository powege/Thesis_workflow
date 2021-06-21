### SCRIPT that calculates the probabilities of SNV given differnet kmer models

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
counts.chr.file <- args[1]
k1CG.file.specific <- args[2]
k3.file.specific <- args[3]
k5.file.specific <- args[4]
k7.file.specific <- args[5]
k9.file.specific <- args[6]
k11.file.specific <- args[7]

# counts.chr.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
# k1CG.file.specific <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz" # outfile
# k3.file.specific <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz" # outfile
# k5.file.specific <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz" # outfile
# k7.file.specific <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz" # outfile
# k9.file.specific <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz" # outfile
# k11.file.specific <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz" # outfile


### IMPORT
k11_specific <- fread(paste0("gunzip -cq ", counts.chr.file))

### FORMAT 

# sum across chromosomes
k11_specific <- setDT(k11_specific)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k11_from, k11_to)] # sum across chromosomes
names(k11_specific)[names(k11_specific) == 'k11_to'] <- 'to'

# 9mer
k9_specific <- k11_specific
colnames(k9_specific) <- gsub("k11", "k9", colnames(k11_specific)) # colnames
k9_specific$k9_from <- stri_sub(k9_specific$k9_from, 2, -2)
k9_specific <- setDT(k9_specific)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, to)] # sum across kmers

# 7mer
k7_specific <- k11_specific
colnames(k7_specific) <- gsub("k11", "k7", colnames(k11_specific)) # colnames
k7_specific$k7_from <- stri_sub(k7_specific$k7_from, 3, -3)
k7_specific <- setDT(k7_specific)[, .(k7_from_N = sum(k7_from_N), k7_to_N = sum(k7_to_N)), by = .(k7_from, to)] # sum across kmers

# 5mer
k5_specific <- k11_specific
colnames(k5_specific) <- gsub("k11", "k5", colnames(k11_specific)) # colnames
k5_specific$k5_from <- stri_sub(k5_specific$k5_from, 4, -4)
k5_specific <- setDT(k5_specific)[, .(k5_from_N = sum(k5_from_N), k5_to_N = sum(k5_to_N)), by = .(k5_from, to)] # sum across kmers

# 3mer
k3_specific <- k11_specific
colnames(k3_specific) <- gsub("k11", "k3", colnames(k11_specific))
k3_specific$k3_from <- stri_sub(k3_specific$k3_from, 5, -5)
k3_specific <- setDT(k3_specific)[, .(k3_from_N = sum(k3_from_N), k3_to_N = sum(k3_to_N)), by = .(k3_from, to)] # sum across kmers

# 1mer
k1_specific <- k11_specific
colnames(k1_specific) <- gsub("k11", "k1", colnames(k11_specific))
k1_specific$k1_from <- stri_sub(k1_specific$k1_from, 6, -6)
k1_specific <- setDT(k1_specific)[, .(k1_from_N = sum(k1_from_N), k1_to_N = sum(k1_to_N)), by = .(k1_from, to)] # sum across kmers

# k1 CG
k1CG_1 <- k3_specific[which(stri_sub(k3_specific$k3_from, 2, -2) == "A" | stri_sub(k3_specific$k3_from, 2, -2) == "C")]
k1CG_1$k1CG_from <- NA
k1CG_1$k1CG_from[ unique(c( grep("CG", k1CG_1$k3_from))) ] <- "C (CG)"
k1CG_1$k1CG_from[is.na(k1CG_1$k1CG_from)] <- "C (nonCG)"
k1CG_1$k1CG_from[which(stri_sub(k1CG_1$k3_from, 2, -2) == "A")] <- "A"
if (length(unique(k1CG_1$k3_from_N)) == 32){
  CG_sub <- k1CG_1[k1CG_from == "C (CG)"]
  nonCG_sub <- k1CG_1[k1CG_from == "C (nonCG)"]
  A_sub <- k1CG_1[k1CG_from == "A"]
  k1CG_2 <- data.table(k1CG_from = c(rep("C (CG)", 3), rep("C (nonCG)", 3), rep("A", 3)),
                       to = c(rep(c("A", "G", "T"), 2), rep(c("C", "G", "T"), 1)),
                       k1CG_from_N = c(rep(sum(unique(CG_sub$k3_from_N)), 3),
                                       rep(sum(unique(nonCG_sub$k3_from_N)), 3),
                                       rep(sum(unique(A_sub$k3_from_N)), 3)),
                       k1CG_to_N =c(sum(CG_sub$k3_to_N[CG_sub$to == "A"]),
                                    sum(CG_sub$k3_to_N[CG_sub$to == "G"]),
                                    sum(CG_sub$k3_to_N[CG_sub$to == "T"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$to == "A"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$to == "G"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$to == "T"]),
                                    sum(A_sub$k3_to_N[A_sub$to == "C"]),
                                    sum(A_sub$k3_to_N[A_sub$to == "G"]),
                                    sum(A_sub$k3_to_N[A_sub$to == "T"]))
  )
  rm(CG_sub, nonCG_sub, A_sub)
}
k1CG_1 <- k1CG_1[,c("k3_from", "to", "k1CG_from")]
k1CG_specific <- k1CG_2[k1CG_1, on = c("k1CG_from", "to")]
rm(k1CG_1, k1CG_2)

# kmer rates of change
k11_specific$k11_mu_rate <- k11_specific$k11_to_N / k11_specific$k11_from_N
k9_specific$k9_mu_rate <- k9_specific$k9_to_N / k9_specific$k9_from_N
k7_specific$k7_mu_rate <- k7_specific$k7_to_N / k7_specific$k7_from_N
k5_specific$k5_mu_rate <- k5_specific$k5_to_N / k5_specific$k5_from_N
k3_specific$k3_mu_rate <- k3_specific$k3_to_N / k3_specific$k3_from_N
k1_specific$k1_mu_rate <- k1_specific$k1_to_N / k1_specific$k1_from_N
k1CG_specific$k1CG_mu_rate <- k1CG_specific$k1CG_to_N / k1CG_specific$k1CG_from_N

k11_specific$k11_mu_rate[k11_specific$k11_from_N == 0] <- 0
k9_specific$k9_mu_rate[k9_specific$k9_from_N == 0] <- 0
k7_specific$k7_mu_rate[k7_specific$k7_from_N == 0] <- 0
k5_specific$k5_mu_rate[k5_specific$k5_from_N == 0] <- 0
k3_specific$k3_mu_rate[k3_specific$k3_from_N == 0] <- 0
k1_specific$k1_mu_rate[k1_specific$k1_from_N == 0] <- 0
k1CG_specific$k1CG_mu_rate[k1CG_specific$k1CG_from_N == 0] <- 0

### OUTPUT
fwrite(k1CG_specific, k1CG.file.specific, compress = "gzip")
fwrite(k3_specific, k3.file.specific, compress = "gzip")
fwrite(k5_specific, k5.file.specific, compress = "gzip")
fwrite(k7_specific, k7.file.specific, compress = "gzip")
fwrite(k9_specific, k9.file.specific, compress = "gzip")
fwrite(k11_specific, k11.file.specific, compress = "gzip")

