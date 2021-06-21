### SCRIPT that calculates the probabilities of SNV for methylated and unmethylated CGs given differnet kmer models

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

#####

# ### RUN ON RESCOMP
# # module load R/3.4.3
# 
# ### SET VARS
# p.mu.file.chr <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_counts_by_chr.csv.gz" # P(mu) 7mer file
# p.mu.file.all <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
# # p.mu.file.chr <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_counts_by_chr.csv.gz" # P(mu) 7mer file
# # p.mu.file.all <- "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
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

### FUNCTIONS 
alakazam <- function(k11_specific, status){
  
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
k11_specific <- k11_specific[k1_from == "C"]
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

# k1CG
k1CG_specific <- k11_specific
colnames(k1CG_specific) <- gsub("k11", "k1CG", colnames(k11_specific))
k1CG_specific$k1_from <- stri_sub(k1CG_specific$k1CG_from, 6, -6)
k1CG_specific <- setDT(k1CG_specific)[, .(k1CG_from_N = sum(k1CG_from_N), k1CG_to_N = sum(k1CG_to_N)), by = .(k1_from, to)] # sum across kmers
k1CG_specific$k1CG_from <- "C (CG)"

# kmer rates of change
k11_specific$k11_mu_rate <- k11_specific$k11_to_N / (k11_specific$k11_from_N)
k9_specific$k9_mu_rate <- k9_specific$k9_to_N / (k9_specific$k9_from_N)
k7_specific$k7_mu_rate <- k7_specific$k7_to_N / (k7_specific$k7_from_N)
k5_specific$k5_mu_rate <- k5_specific$k5_to_N / (k5_specific$k5_from_N)
k3_specific$k3_mu_rate <- k3_specific$k3_to_N / (k3_specific$k3_from_N)
k1CG_specific$k1CG_mu_rate <- k1CG_specific$k1CG_to_N / (k1CG_specific$k1CG_from_N)

k11_specific$k11_mu_rate[k11_specific$k11_from_N == 0] <- 0
k9_specific$k9_mu_rate[k9_specific$k9_from_N == 0] <- 0
k7_specific$k7_mu_rate[k7_specific$k7_from_N == 0] <- 0
k5_specific$k5_mu_rate[k5_specific$k5_from_N == 0] <- 0
k3_specific$k3_mu_rate[k3_specific$k3_from_N == 0] <- 0
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
output_specific <- k1CG_specific[output_specific, on = c("k1_from", "to")]

output_specific$meth_status <- status

# sum by k_from
k3_any <- setDT(k3_specific)[, .(k3_to_N = sum(k3_to_N)), by = .(k3_from, k3_from_N)]
k5_any <- setDT(k5_specific)[, .(k5_to_N = sum(k5_to_N)), by = .(k5_from, k5_from_N)]
k7_any <- setDT(k7_specific)[, .(k7_to_N = sum(k7_to_N)), by = .(k7_from, k7_from_N)]
k9_any <- setDT(k9_specific)[, .(k9_to_N = sum(k9_to_N)), by = .(k9_from, k9_from_N)]
k11_any <- setDT(k11_specific)[, .(k11_to_N = sum(k11_to_N)), by = .(k11_from, k11_from_N)]
k1CG_any <- setDT(k1CG_specific)[, .(k1CG_to_N = sum(k1CG_to_N)), by = .(k1CG_from, k1_from, k1CG_from_N)]

# kmer rates of change
k11_any$k11_mu_rate <- k11_any$k11_to_N / k11_any$k11_from_N
k9_any$k9_mu_rate <- k9_any$k9_to_N / k9_any$k9_from_N
k7_any$k7_mu_rate <- k7_any$k7_to_N / k7_any$k7_from_N
k5_any$k5_mu_rate <- k5_any$k5_to_N / k5_any$k5_from_N
k3_any$k3_mu_rate <- k3_any$k3_to_N / k3_any$k3_from_N
k1CG_any$k1CG_mu_rate <- k1CG_any$k1CG_to_N / k1CG_any$k1CG_from_N

k11_any$k11_mu_rate[k11_any$k11_from_N == 0] <- 0
k9_any$k9_mu_rate[k9_any$k9_from_N == 0] <- 0
k7_any$k7_mu_rate[k7_any$k7_from_N == 0] <- 0
k5_any$k5_mu_rate[k5_any$k5_from_N == 0] <- 0
k3_any$k3_mu_rate[k3_any$k3_from_N == 0] <- 0
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
output_any <- k1CG_any[output_any, on = c("k1_from")]

output_any$meth_status <- status

return(list(output_specific, output_any))
}


### SET VARS
p.mu.file.all.meth <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
p.mu.file.all.unmeth <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
out.file.specific.meth <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_methylated_kmer_pSNV_specific.csv.gz" # outfile
out.file.specific.unmeth <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_kmer_pSNV_specific.csv.gz" # outfile
out.file.any.meth <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_methylated_kmer_pSNV_any.csv.gz" # outfile
out.file.any.unmeth <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_kmer_pSNV_any.csv.gz" # outfile

### IMPORT
k11_meth <- fread(paste0("gunzip -cq ", p.mu.file.all.meth))
k11_unmeth <- fread(paste0("gunzip -cq ", p.mu.file.all.unmeth))

### RUN
meth_list <- alakazam(k11_specific = k11_meth, status = "methylated")
unmeth_list <- alakazam(k11_specific = k11_unmeth, status = "unmethylated")

### FORMAT
out_specific_meth <- meth_list[[1]]
out_specific_unmeth <- unmeth_list[[1]]

### EXPORT
fwrite(out_specific_unmeth, out.file.specific.unmeth, compress = "gzip")
fwrite(out_specific_meth, out.file.specific.meth, compress = "gzip")

#####
