rm(list = ls())
graphics.off()

library(data.table)

### IMPORT
k11_specific <- fread(paste0("gunzip -cq ", "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_autosome_counts.csv.gz"))


### FORMAT

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

### FUNCTION
alakazam <- function(k11_specific){
  
colnames(k11_specific) <- c("k11_from", "k11_to", "k11_from_N", "k11_to_N")
  
# 11mer
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

# kmer rates of change
k11_specific$k11_to_rate <- k11_specific$k11_to_N / k11_specific$k11_from_N
k9_specific$k9_to_rate <- k9_specific$k9_to_N / k9_specific$k9_from_N
k7_specific$k7_to_rate <- k7_specific$k7_to_N / k7_specific$k7_from_N
k5_specific$k5_to_rate <- k5_specific$k5_to_N / k5_specific$k5_from_N
k3_specific$k3_to_rate <- k3_specific$k3_to_N / k3_specific$k3_from_N
k1_specific$k1_to_rate <- k1_specific$k1_to_N / k1_specific$k1_from_N

k11_specific$k11_to_rate[k11_specific$k11_from_N == 0] <- 0
k9_specific$k9_to_rate[k9_specific$k9_from_N == 0] <- 0
k7_specific$k7_to_rate[k7_specific$k7_from_N == 0] <- 0
k5_specific$k5_to_rate[k5_specific$k5_from_N == 0] <- 0
k3_specific$k3_to_rate[k3_specific$k3_from_N == 0] <- 0
k1_specific$k1_to_rate[k1_specific$k1_from_N == 0] <- 0

out <- list(k11_specific,
            k9_specific,
            k7_specific,
            k5_specific,
            k3_specific,
            k1_specific)
return(out)
}

mod_5 <- alakazam(k11_specific)
mod_3 <- alakazam(k11_complemnt)

kmers <- c("k11", "k9", "k7", "k5", "k3", "k1")
estimate <- rep(NA, 6)
df <- rep(NA, 6)
p_val <- rep(NA, 6)
method <- rep(NA, 6)
for (i in 1:length(kmers)){
  
  dt5 <- mod_5[[i]]
  dt3 <- mod_3[[i]]
  
  colnames(dt5) <- c("from", "to", "from_N_5", "to_N_5", "mu_rate_5")
  colnames(dt3) <- c("from", "to", "from_N_3", "to_N_3", "mu_rate_3")
  dt <- dt5[dt3, on = c("from", "to")]
  cor <- cor.test(dt$mu_rate_5, dt$mu_rate_3, method = "pearson")
  estimate[i] <- cor$estimate
  df[i] <- cor$parameter
  p_val[i] <- cor$p.val
  method[i] <- cor$method
  print(i)
}

output <- data.table(kmers, 
                     estimate,
                     df,
                     p_val,
                     method)
### OUTPUT
fwrite(output, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Figures_tables/Table_complementary_kmer_cor.csv")

