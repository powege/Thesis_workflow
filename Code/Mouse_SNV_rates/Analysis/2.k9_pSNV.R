### SCRIPT that calculates the probabilities of SNV given differnet kmer models

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)
library(plyr)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 

### SET VARS
counts.chr.file <- args[1]
k1.outfile <- args[2]
k1CG.outfile <- args[3]
k3.outfile <- args[4]
k5.outfile <- args[5]
k7.outfile <- args[6]
k9.outfile <- args[7]

# counts.chr.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k9_SNV_counts_by_chr.csv.gz" # P(mu)
# k1.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1_pSNV.csv.gz" # outfile
# k1CG.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_pSNV.csv.gz" # outfile
# k3.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_pSNV.csv.gz" # outfile
# k5.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k5_pSNV.csv.gz" # outfile
# k7.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV.csv.gz" # outfile
# k9.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k9_pSNV.csv.gz" # outfile

### FUNCTIONS

# FUNCTION that sums across complimantary strands (ie AGC == GCT)
complement <- function(forward){
  
  colnames(forward) <- c("k9_from", "k9_to", "k9_from_N", "k9_to_N")
  complement <- forward # reverse strand
  complement$k9_from <- stri_reverse(complement$k9_from)
  complement$k9_to <- stri_reverse(complement$k9_to)
  
  # replace all bases with complement
  complement$k9_from <- gsub("A", "B", complement$k9_from)
  complement$k9_from <- gsub("C", "D", complement$k9_from)
  complement$k9_from <- gsub("T", "A", complement$k9_from)
  complement$k9_from <- gsub("G", "C", complement$k9_from)
  complement$k9_from <- gsub("B", "T", complement$k9_from)
  complement$k9_from <- gsub("D", "G", complement$k9_from)
  complement$k9_to <- gsub("A", "B", complement$k9_to)
  complement$k9_to <- gsub("C", "D", complement$k9_to)
  complement$k9_to <- gsub("T", "A", complement$k9_to)
  complement$k9_to <- gsub("G", "C", complement$k9_to)
  complement$k9_to <- gsub("B", "T", complement$k9_to)
  complement$k9_to <- gsub("D", "G", complement$k9_to)
  
  # sum across complements
  output <- rbind(forward, complement)
  output <- setDT(output)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, k9_to)]
  
  return(output)
}

### IMPORT
k9_specific <- fread(paste0("gunzip -cq ", counts.chr.file))

### FORMAT 
colnames(k9_specific) <- c("chromosome", "k9_from", "k9_to", "k9_from_N", "k9_to_N")

# sum across chromosomes
k9_specific <- setDT(k9_specific)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, k9_to)] # sum across chromosomes

# sum across complementary kmers
k9_specific <- complement(k9_specific)

# calculate k1 change
k9_specific$k1_from <- stri_sub(k9_specific$k9_from, 5, -5)
k9_specific$k1_to <- stri_sub(k9_specific$k9_to, 5, -5)

# subset A and C mutations
k9_specific <- k9_specific[k1_from %in% c("A", "C")]

# 7mer
k7_specific <- k9_specific
colnames(k7_specific) <- gsub("k9", "k7", colnames(k9_specific)) # colnames
k7_specific$k7_from <- stri_sub(k7_specific$k7_from, 2, -2)
k7_specific$k7_to <- stri_sub(k7_specific$k7_to, 2, -2)
k7_specific <- setDT(k7_specific)[, .(k7_from_N = sum(k7_from_N), k7_to_N = sum(k7_to_N)), by = .(k7_from, k7_to, k1_from, k1_to)] # sum across kmers

# 5mer
k5_specific <- k9_specific
colnames(k5_specific) <- gsub("k9", "k5", colnames(k9_specific)) # colnames
k5_specific$k5_from <- stri_sub(k5_specific$k5_from, 3, -3)
k5_specific$k5_to <- stri_sub(k5_specific$k5_to, 3, -3)
k5_specific <- setDT(k5_specific)[, .(k5_from_N = sum(k5_from_N), k5_to_N = sum(k5_to_N)), by = .(k5_from, k5_to, k1_from, k1_to)] # sum across kmers

# 3mer
k3_specific <- k9_specific
colnames(k3_specific) <- gsub("k9", "k3", colnames(k9_specific))
k3_specific$k3_from <- stri_sub(k3_specific$k3_from, 4, -4)
k3_specific$k3_to <- stri_sub(k3_specific$k3_to, 4, -4)
k3_specific <- setDT(k3_specific)[, .(k3_from_N = sum(k3_from_N), k3_to_N = sum(k3_to_N)), by = .(k3_from, k3_to, k1_from, k1_to)] # sum across kmers

# 1mer
k1_specific <- k9_specific[,c("k1_from", "k1_to", "k9_from_N", "k9_to_N")]
colnames(k1_specific) <- gsub("k9", "k1", colnames(k1_specific))
k1_specific <- setDT(k1_specific)[, .(k1_from_N = sum(k1_from_N), k1_to_N = sum(k1_to_N)), by = .(k1_from, k1_to)] # sum across kmers

# k1 CG
k1CG_1 <- k3_specific[k3_specific$k1_from == "A" | k3_specific$k1_from == "C"]
k1CG_1$k1CG_group <- NA
k1CG_1$k1CG_group[ unique(c( grep("CG", k1CG_1$k3_from))) ] <- "C (CG)"
k1CG_1$k1CG_group[is.na(k1CG_1$k1CG_group)] <- "C (nonCG)"
k1CG_1$k1CG_group[which(stri_sub(k1CG_1$k3_from, 2, -2) == "A")] <- "A"
if (length(unique(k1CG_1$k3_from_N)) == 32){
  CG_sub <- k1CG_1[k1CG_group == "C (CG)"]
  nonCG_sub <- k1CG_1[k1CG_group == "C (nonCG)"]
  A_sub <- k1CG_1[k1CG_group == "A"]
  k1CG_2 <- data.table(k1CG_group = c(rep("C (CG)", 3), rep("C (nonCG)", 3), rep("A", 3)),
                       k1_to = c(rep(c("A", "G", "T"), 2), rep(c("C", "G", "T"), 1)),
                       k1CG_from_N = c(rep(sum(unique(CG_sub$k3_from_N)), 3),
                                       rep(sum(unique(nonCG_sub$k3_from_N)), 3),
                                       rep(sum(unique(A_sub$k3_from_N)), 3)),
                       k1CG_to_N =c(sum(CG_sub$k3_to_N[CG_sub$k1_to == "A"]),
                                    sum(CG_sub$k3_to_N[CG_sub$k1_to == "G"]),
                                    sum(CG_sub$k3_to_N[CG_sub$k1_to == "T"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$k1_to == "A"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$k1_to == "G"]),
                                    sum(nonCG_sub$k3_to_N[nonCG_sub$k1_to == "T"]),
                                    sum(A_sub$k3_to_N[A_sub$k1_to == "C"]),
                                    sum(A_sub$k3_to_N[A_sub$k1_to == "G"]),
                                    sum(A_sub$k3_to_N[A_sub$k1_to == "T"]))
  )
  rm(CG_sub, nonCG_sub, A_sub)
}
k1CG_1 <- k1CG_1[,c("k3_from", "k3_to", "k1_from", "k1_to", "k1CG_group")]
k1CG_specific <- k1CG_2[k1CG_1, on = c("k1CG_group", "k1_to")]
rm(k1CG_1, k1CG_2)

# any substitution
k9_any <- setDT(k9_specific)[, .(k9_to_N = sum(k9_to_N)), by = .(k9_from, k1_from, k9_from_N)] # sum across kmers
k7_any <- setDT(k7_specific)[, .(k7_to_N = sum(k7_to_N)), by = .(k7_from, k1_from, k7_from_N)] # sum across kmers
k5_any <- setDT(k5_specific)[, .(k5_to_N = sum(k5_to_N)), by = .(k5_from, k1_from, k5_from_N)] # sum across kmers
k3_any <- setDT(k3_specific)[, .(k3_to_N = sum(k3_to_N)), by = .(k3_from, k1_from, k3_from_N)] # sum across kmers
k1_any <- setDT(k1_specific)[, .(k1_to_N = sum(k1_to_N)), by = .(k1_from, k1_from_N)] # sum across kmers
k1CG_any <- setDT(k1CG_specific)[, .(k1CG_to_N = sum(k1CG_to_N)), by = .(k1CG_group, k1_from, k3_from, k1CG_from_N)] # sum across kmers

k9_any$k1_to <- "*"
k7_any$k1_to <- "*"
k5_any$k1_to <- "*"
k3_any$k1_to <- "*"
k1_any$k1_to <- "*"
k1CG_any$k1_to <- "*"

k9_out <- rbind.fill(k9_specific, k9_any)
k7_out <- rbind.fill(k7_specific, k7_any)
k5_out <- rbind.fill(k5_specific, k5_any)
k3_out <- rbind.fill(k3_specific, k3_any)
k1_out <- rbind.fill(k1_specific, k1_any)
k1CG_out <- rbind.fill(k1CG_specific, k1CG_any)

# kmer rates of change
k9_out$k9_mu_rate <- k9_out$k9_to_N / k9_out$k9_from_N
k7_out$k7_mu_rate <- k7_out$k7_to_N / k7_out$k7_from_N
k5_out$k5_mu_rate <- k5_out$k5_to_N / k5_out$k5_from_N
k3_out$k3_mu_rate <- k3_out$k3_to_N / k3_out$k3_from_N
k1_out$k1_mu_rate <- k1_out$k1_to_N / k1_out$k1_from_N
k1CG_out$k1CG_mu_rate <- k1CG_out$k1CG_to_N / k1CG_out$k1CG_from_N

k9_out$k9_mu_rate[k9_out$k9_from_N == 0] <- 0
k9_out$k9_mu_rate[k9_out$k9_from_N == 0] <- 0
k7_out$k7_mu_rate[k7_out$k7_from_N == 0] <- 0
k5_out$k5_mu_rate[k5_out$k5_from_N == 0] <- 0
k3_out$k3_mu_rate[k3_out$k3_from_N == 0] <- 0
k1_out$k1_mu_rate[k1_out$k1_from_N == 0] <- 0
k1CG_out$k1CG_mu_rate[k1CG_out$k1CG_from_N == 0] <- 0

### OUTPUT
fwrite(k1CG_out, k1CG.outfile, compress = "gzip")
fwrite(k1_out, k1.outfile, compress = "gzip")
fwrite(k3_out, k3.outfile, compress = "gzip")
fwrite(k5_out, k5.outfile, compress = "gzip")
fwrite(k7_out, k7.outfile, compress = "gzip")
fwrite(k9_out, k9.outfile, compress = "gzip")

#####

# hist(k3_specific$k3_mu_rate[which(stri_sub(k3_specific$k3_from, 2, -1) == "CG" & k3_specific$k1_to == "T")])
# hist(k5_specific$k5_mu_rate[which(stri_sub(k5_specific$k5_from, 3, -2) == "CG" & k5_specific$k1_to == "T")])
# hist(k7_specific$k7_mu_rate[which(stri_sub(k7_specific$k7_from, 4, -3) == "CG" & k7_specific$k1_to == "T")])
# hist(k9_specific$k9_mu_rate[which(stri_sub(k9_specific$k9_from, 5, -4) == "CG" & k9_specific$k1_to == "T")])















