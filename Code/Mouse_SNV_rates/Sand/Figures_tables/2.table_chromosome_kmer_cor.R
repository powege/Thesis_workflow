rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

### IMPORT 
k11_specific <- fread(paste0("gunzip -cq ", "/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_counts_by_chr.csv.gz"))

### FORMAT
colnames(k11_specific) <- c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome") # colnames

### SUM COMPLEMENTS 
k11_complemnt <- k11_specific # reverse strand
k11_complemnt$k11_from <- stri_reverse(k11_complemnt$k11_from)
k11_complemnt$k11_to <- stri_reverse(k11_complemnt$k11_to)

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
k11_specific <- setDT(k11_specific)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k11_from, k11_to, chromosome)]
rm(k11_complemnt)

### SUBSET A & C
k11_specific$k1_from <- stri_sub(k11_specific$k11_from, 6, -6)
k11_specific <- k11_specific[k1_from == "A" | k1_from == "C"]
k11_specific$to <- stri_sub(k11_specific$k11_to, 6, -6)
k11_specific <- k11_specific[,c("k11_from", "to", "k11_from_N", "k11_to_N", "chromosome")]

k11_list <- list()
k9_list <- list()
k7_list <- list()
k5_list <- list()
k3_list <- list()
k1_list <- list()
for (chr in 1:19){
  
  k11_chr <- k11_specific[chromosome == chr]
  k11_chr$chromosome <- NULL
  
  # 9mer
  k9_chr <- k11_chr
  colnames(k9_chr) <- gsub("k11", "k9", colnames(k11_chr)) # colnames
  k9_chr$k9_from <- stri_sub(k9_chr$k9_from, 2, -2)
  k9_chr <- setDT(k9_chr)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, to)] # sum across kmers
  
  # 7mer
  k7_chr <- k11_chr
  colnames(k7_chr) <- gsub("k11", "k7", colnames(k11_chr)) # colnames
  k7_chr$k7_from <- stri_sub(k7_chr$k7_from, 3, -3)
  k7_chr <- setDT(k7_chr)[, .(k7_from_N = sum(k7_from_N), k7_to_N = sum(k7_to_N)), by = .(k7_from, to)] # sum across kmers
  
  # 5mer
  k5_chr <- k11_chr
  colnames(k5_chr) <- gsub("k11", "k5", colnames(k11_chr)) # colnames
  k5_chr$k5_from <- stri_sub(k5_chr$k5_from, 4, -4)
  k5_chr <- setDT(k5_chr)[, .(k5_from_N = sum(k5_from_N), k5_to_N = sum(k5_to_N)), by = .(k5_from, to)] # sum across kmers
  
  # 3mer
  k3_chr <- k11_chr
  colnames(k3_chr) <- gsub("k11", "k3", colnames(k11_chr))
  k3_chr$k3_from <- stri_sub(k3_chr$k3_from, 5, -5)
  k3_chr <- setDT(k3_chr)[, .(k3_from_N = sum(k3_from_N), k3_to_N = sum(k3_to_N)), by = .(k3_from, to)] # sum across kmers
  
  # 1mer
  k1_chr <- k11_chr
  colnames(k1_chr) <- gsub("k11", "k1", colnames(k11_chr))
  k1_chr$k1_from <- stri_sub(k1_chr$k1_from, 6, -6)
  k1_chr <- setDT(k1_chr)[, .(k1_from_N = sum(k1_from_N), k1_to_N = sum(k1_to_N)), by = .(k1_from, to)] # sum across kmers
  
  # kmer rates of change
  k11_to_rate <- k11_chr$k11_to_N / k11_chr$k11_from_N
  k9_to_rate <- k9_chr$k9_to_N / k9_chr$k9_from_N
  k7_to_rate <- k7_chr$k7_to_N / k7_chr$k7_from_N
  k5_to_rate <- k5_chr$k5_to_N / k5_chr$k5_from_N
  k3_to_rate <- k3_chr$k3_to_N / k3_chr$k3_from_N
  k1_to_rate <- k1_chr$k1_to_N / k1_chr$k1_from_N

  k11_to_rate[k11_chr$k11_from_N == 0] <- 0
  k9_to_rate[k9_chr$k9_from_N == 0] <- 0
  k7_to_rate[k7_chr$k7_from_N == 0] <- 0
  k5_to_rate[k5_chr$k5_from_N == 0] <- 0
  k3_to_rate[k3_chr$k3_from_N == 0] <- 0
  k1_to_rate[k1_chr$k1_from_N == 0] <- 0
  
  k11_list[[chr]] <- k11_to_rate
  k9_list[[chr]] <- k9_to_rate
  k7_list[[chr]] <- k7_to_rate
  k5_list[[chr]] <- k5_to_rate
  k3_list[[chr]] <- k3_to_rate
  k1_list[[chr]] <- k1_to_rate
  
  print(chr)
}

k11_dt <- do.call("cbind", k11_list)
k9_dt <- do.call("cbind", k9_list)
k7_dt <- do.call("cbind", k7_list)
k5_dt <- do.call("cbind", k5_list)
k3_dt <- do.call("cbind", k3_list)
k1_dt <- do.call("cbind", k1_list)

k11_out <- cor(k11_dt)
k9_out <- cor(k9_dt)
k7_out <- cor(k7_dt)
k5_out <- cor(k5_dt)
k3_out <- cor(k3_dt)
k1_out <- cor(k1_dt)

k11_out <- as.data.table(k11_out)
k9_out <- as.data.table(k9_out)
k7_out <- as.data.table(k7_out)
k5_out <- as.data.table(k5_out)
k3_out <- as.data.table(k3_out)
k1_out <- as.data.table(k1_out)

### EXPORT
fwrite(k11_out, "/well/lindgren/George/Data/Thesis_workflow/Results/Figures_tables/table_k11_pSNV_chromosome_cor.csv")
fwrite(k9_out, "/well/lindgren/George/Data/Thesis_workflow/Results/Figures_tables/table_k9_pSNV_chromosome_cor.csv")
fwrite(k7_out, "/well/lindgren/George/Data/Thesis_workflow/Results/Figures_tables/table_k7_pSNV_chromosome_cor.csv")
fwrite(k5_out, "/well/lindgren/George/Data/Thesis_workflow/Results/Figures_tables/table_k5_pSNV_chromosome_cor.csv")
fwrite(k3_out, "/well/lindgren/George/Data/Thesis_workflow/Results/Figures_tables/table_k3_pSNV_chromosome_cor.csv")
fwrite(k1_out, "/well/lindgren/George/Data/Thesis_workflow/Results/Figures_tables/table_k1_pSNV_chromosome_cor.csv")

