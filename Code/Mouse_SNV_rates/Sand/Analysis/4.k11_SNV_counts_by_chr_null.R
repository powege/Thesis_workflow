### SCRIPT that creates null datasets with same k11_from_N as observed:
# 1. where k11_to_N is sampled from a binomial distribution with same pMU as k1CG model
# ... k3 - k9
# 2. where k11_to_N is sampled from a binomial distribution with same pMU as k11 model 

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

# FUNCTION that sums across complimantary strands (ie AGC == GCT)
complement <- function(forward){
  
  complement <- forward # reverse strand
  complement$k11_from <- stri_reverse(complement$k11_from)
  
  # replace all bases with complement
  complement$k11_from <- gsub("A", "B", complement$k11_from)
  complement$k11_from <- gsub("C", "D", complement$k11_from)
  complement$k11_from <- gsub("T", "A", complement$k11_from)
  complement$k11_from <- gsub("G", "C", complement$k11_from)
  complement$k11_from <- gsub("B", "T", complement$k11_from)
  complement$k11_from <- gsub("D", "G", complement$k11_from)
  complement$to <- gsub("A", "B", complement$to)
  complement$to <- gsub("C", "D", complement$to)
  complement$to <- gsub("T", "A", complement$to)
  complement$to <- gsub("G", "C", complement$to)
  complement$to <- gsub("B", "T", complement$to)
  complement$to <- gsub("D", "G", complement$to)
  
  # sum across complements
  output <- rbind(forward, complement)

  return(output)
}

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 

### SET VARS
p.mu.chr.file <- args[1]
k11.p.snv.specific.file <- args[2]
k11.null.outfile <- args[3]
k1CG.null.outfile <- args[4]

# p.mu.chr.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
# k11.p.snv.specific.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_methylated_kmer_pSNV_specific.csv.gz"
# k11.null.outfile <- ""
# k1CG.null.outfile <- ""

# p.mu.chr.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
# k11.p.snv.specific.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz"
# k11.null.outfile <- ""
# k1CG.null.outfile <- ""

### IMPORT 
p_snv_specific <- fread(paste0("gunzip -cq ", k11.p.snv.specific.file))
k11_chr <- fread(paste0("gunzip -cq ", p.mu.chr.file))
# tmp <- k11_chr
# tmp$chromosome <- 18
# k11_chr$chromosome <- 19
# k11_chr <- rbind(k11_chr, tmp)
# k11_chr <- k11_chr[,c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome")] # colnames
# rm(tmp)

### FORMAT
colnames(k11_chr) <- c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome") # colnames
k11_chr$to <- stri_sub(k11_chr$k11_to, 6, -6)
k11_chr <- k11_chr[,c("k11_from", "k11_to", "to", "k11_from_N", "chromosome")]

k1CG <- unique(p_snv_specific[,c("k11_from", "to", "k1CG_mu_rate")])
k11 <- unique(p_snv_specific[,c("k11_from", "to", "k11_mu_rate")])
k11$k11_mu_rate[k11$k11_mu_rate == 0] <- min(k11$k11_mu_rate[k11$k11_mu_rate > 0])
rm(p_snv_specific)

k1CG <- complement(k1CG)
k11 <- complement(k11)

k1CG_null <- k11_chr[k1CG, on = c("k11_from", "to")]
k11_null <- k11_chr[k11, on = c("k11_from", "to")]
rm(k1CG, k11, k11_chr)

k1CG_null$k11_to_N <- rbinom(size = k1CG_null$k11_from_N, n = nrow(k1CG_null), prob = k1CG_null$k1CG_mu_rate)
k11_null$k11_to_N <- rbinom(size = k11_null$k11_from_N, n = nrow(k11_null), prob = k11_null$k11_mu_rate)

k11_null <- k11_null[,c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome")]
k1CG_null <- k1CG_null[,c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome")]

### EXPORT 
fwrite(k11_null, k11.null.outfile, compress = "gzip")
fwrite(k1CG_null, k1CG.null.outfile, compress = "gzip")






