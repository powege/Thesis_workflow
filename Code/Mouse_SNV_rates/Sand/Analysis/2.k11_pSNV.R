### SCRIPT that calculates the probabilities of SNV given differnet kmer models

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

#####

# ### RUN ON RESCOMP
# # module load R/3.4.3
# 
# ### SET VARS
# p.mu.file.chr <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_counts_by_chr.csv.gz" # P(mu) 7mer file
# p.mu.file.all <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
# 
# ### IMPORT
# k11_chr <- fread(paste0("gunzip -cq ", p.mu.file.chr))
# 
# ### FORMAT
# colnames(k11_chr) <- c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome") # colnames
# k11_specific <- setDT(k11_chr)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k11_from, k11_to)] # sum across chromosomes
# 
# ### EXPORT
# fwrite(k11_specific, p.mu.file.all, compress = "gzip")

#####

### SET VARS
p.mu.file.all <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
out.file.specific <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz" # outfile
out.file.any <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_any.csv.gz" # outfile

### IMPORT
k11_specific <- fread(paste0("gunzip -cq ", p.mu.file.all))

## average complementary strands

# reverse strand
k11_complemnt <- k11_specific
k11_complemnt$k11_from <- stri_reverse(k11_complemnt$k11_from)
k11_complemnt$k11_to <- stri_reverse(k11_complemnt$k11_to)

# replace all bases with complement
k11_complemnt$k11_from <- gsub("A", "B", k11_complemnt$k11_from)
k11_complemnt$k11_from <- gsub("C", "D", k11_complemnt$k11_from)
k11_complemnt$k11_from <- gsub("T", "A", k11_complemnt$k11_from)
k11_complemnt$k11_from <- gsub("G", "C", k11_complemnt$k11_from)
k11_complemnt$k11_from <- gsub("B", "T", k11_complemnt$k11_from)
k11_complemnt$k11_from <- gsub("D", "G", k11_complemnt$k11_from)
k11_complemnt$k11_to <- gsub("A", "B", k11_complemnt$k11_to)
k11_complemnt$k11_to <- gsub("C", "D", k11_complemnt$k11_to)
k11_complemnt$k11_to <- gsub("T", "A", k11_complemnt$k11_to)
k11_complemnt$k11_to <- gsub("G", "C", k11_complemnt$k11_to)
k11_complemnt$k11_to <- gsub("B", "T", k11_complemnt$k11_to)
k11_complemnt$k11_to <- gsub("D", "G", k11_complemnt$k11_to)

# sum across complements
k11_specific <- rbind(k11_specific, k11_complemnt)
k11_specific <- setDT(k11_specific)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k11_from, k11_to)]
rm(k11_complemnt)

# subset A and C mutations
k11_specific$k1_from <- stri_sub(k11_specific$k11_from, 6, -6)
k11_specific <- k11_specific[k1_from == "A" | k1_from == "C"]
k11_specific$to <- stri_sub(k11_specific$k11_to, 6, -6)
k11_specific <- k11_specific[,c("k11_from", "to", "k11_from_N", "k11_to_N")]

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

# merge
output_specific <- k11_specific
output_specific$k9_from <- stri_sub(output_specific$k11_from, 2, -2)
output_specific$k7_from <- stri_sub(output_specific$k11_from, 3, -3)
output_specific$k5_from <- stri_sub(output_specific$k11_from, 4, -4)
output_specific$k3_from <- stri_sub(output_specific$k11_from, 5, -5)
output_specific$k1_from <- stri_sub(output_specific$k11_from, 6, -6)
output_specific <- k9_specific[output_specific, on = c("k9_from", "to")]
output_specific <- k7_specific[output_specific, on = c("k7_from", "to")]
output_specific <- k5_specific[output_specific, on = c("k5_from", "to")]
output_specific <- k3_specific[output_specific, on = c("k3_from", "to")]
output_specific <- k1_specific[output_specific, on = c("k1_from", "to")]
output_specific <- k1CG_specific[output_specific, on = c("k3_from", "to")]

# sum by k_from
k1_any <- setDT(k1_specific)[, .(k1_to_N = sum(k1_to_N)), by = .(k1_from, k1_from_N)]
k3_any <- setDT(k3_specific)[, .(k3_to_N = sum(k3_to_N)), by = .(k3_from, k3_from_N)]
k5_any <- setDT(k5_specific)[, .(k5_to_N = sum(k5_to_N)), by = .(k5_from, k5_from_N)]
k7_any <- setDT(k7_specific)[, .(k7_to_N = sum(k7_to_N)), by = .(k7_from, k7_from_N)]
k9_any <- setDT(k9_specific)[, .(k9_to_N = sum(k9_to_N)), by = .(k9_from, k9_from_N)]
k11_any <- setDT(k11_specific)[, .(k11_to_N = sum(k11_to_N)), by = .(k11_from, k11_from_N)]
k1CG_any <- setDT(k1CG_specific)[, .(k1CG_to_N = sum(k1CG_to_N)), by = .(k3_from, k1CG_from_N)]

# kmer rates of change
k11_any$k11_mu_rate <- k11_any$k11_to_N / k11_any$k11_from_N
k9_any$k9_mu_rate <- k9_any$k9_to_N / k9_any$k9_from_N
k7_any$k7_mu_rate <- k7_any$k7_to_N / k7_any$k7_from_N
k5_any$k5_mu_rate <- k5_any$k5_to_N / k5_any$k5_from_N
k3_any$k3_mu_rate <- k3_any$k3_to_N / k3_any$k3_from_N
k1_any$k1_mu_rate <- k1_any$k1_to_N / k1_any$k1_from_N
k1CG_any$k1CG_mu_rate <- k1CG_any$k1CG_to_N / k1CG_any$k1CG_from_Ns

k11_any$k11_mu_rate[k11_any$k11_from_N == 0] <- 0
k9_any$k9_mu_rate[k9_any$k9_from_N == 0] <- 0
k7_any$k7_mu_rate[k7_any$k7_from_N == 0] <- 0
k5_any$k5_mu_rate[k5_any$k5_from_N == 0] <- 0
k3_any$k3_mu_rate[k3_any$k3_from_N == 0] <- 0
k1_any$k1_mu_rate[k1_any$k1_from_N == 0] <- 0
k1CG_any$k1CG_mu_rate[k1CG_any$k1CG_from_N == 0] <- 0

# merge
output_any <- k11_any
output_any$k9_from <- stri_sub(output_any$k11_from, 2, -2)
output_any$k7_from <- stri_sub(output_any$k11_from, 3, -3)
output_any$k5_from <- stri_sub(output_any$k11_from, 4, -4)
output_any$k3_from <- stri_sub(output_any$k11_from, 5, -5)
output_any$k1_from <- stri_sub(output_any$k11_from, 6, -6)
output_any <- k9_any[output_any, on = c("k9_from")]
output_any <- k7_any[output_any, on = c("k7_from")]
output_any <- k5_any[output_any, on = c("k5_from")]
output_any <- k3_any[output_any, on = c("k3_from")]
output_any <- k1_any[output_any, on = c("k1_from")]
output_any <- k1CG_any[output_any, on = c("k3_from")]

### OUTPUT
fwrite(output_specific, out.file.specific, compress = "gzip")
fwrite(output_any, out.file.any, compress = "gzip")

#####

# x <- CG_specific[k3_specific, on = c("k3_from", "to")]
# x$k1_from <- stri_sub(x$k3_from, 2, -2)
# x <- k1_specific[x, on = c("to", "k1_from")]
# x <- x[,c("k1_from", "k3_from", "CG", "to", "k1_mu_rate", "k3_mu_rate", "CG_mu_rate")]



# k11_chr <- fread("/well/lindgren/George/Data/Mu_rates/MGP_v5_allSTRAIN_k11_SNV_counts_by_chr.csv")
# colnames(k11_chr) <- c("k11_from", "k11_from_N", "k11_to", "k11_mu_N", "chromosome") # colnames
# k11_specific <- setDT(k11_chr)[, .(k11_from_N = sum(k11_from_N), k11_mu_N = sum(k11_mu_N)), by = .(k11_from, k11_to)] # sum across chromosomes
# fwrite(k11_specific, "/well/lindgren/George/Data/Mu_rates/MGP_v5_allSTRAIN_k11_SNV_counts_all_chr.csv")


# x <- k11_specific$k11_mu_N / k11_specific$k11_from_N
# y <- rbinom(n = 1:nrow(k11_specific), size = k11_specific$k11_from_N, prob = mean(x, na.rm = T))
# range(y)
# length(y[which(y == 0)])
