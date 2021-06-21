### SCRIPT that runs cross validation by chromosome (selects 4 random chromosomes for the test dataset)

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

### FUNCTIONS 

# FUNCTION that sums across complimantary strands (ie AGC == GCT)
complement <- function(forward){
  
  colnames(forward) <- c("k_from", "k_to", "k_from_N", "k_to_N")
  complement <- forward # reverse strand
  complement$k_from <- stri_reverse(complement$k_from)
  
  # replace all bases with complement
  complement$k_from <- gsub("A", "B", complement$k_from)
  complement$k_from <- gsub("C", "D", complement$k_from)
  complement$k_from <- gsub("T", "A", complement$k_from)
  complement$k_from <- gsub("G", "C", complement$k_from)
  complement$k_from <- gsub("B", "T", complement$k_from)
  complement$k_from <- gsub("D", "G", complement$k_from)
  complement$k_to <- gsub("A", "B", complement$k_to)
  complement$k_to <- gsub("C", "D", complement$k_to)
  complement$k_to <- gsub("T", "A", complement$k_to)
  complement$k_to <- gsub("G", "C", complement$k_to)
  complement$k_to <- gsub("B", "T", complement$k_to)
  complement$k_to <- gsub("D", "G", complement$k_to)
  
  # sum across complements
  output <- rbind(forward, complement)
  output <- setDT(output)[, .(k_from_N = sum(k_from_N), k_to_N = sum(k_to_N)), by = .(k_from, k_to)]
  
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
p.mu.chr.methylated.file <- args[2]
p.mu.chr.unmethylated.file <- args[3]
out.specific.meth.file <- args[4]
out.specific.unmeth.file <- args[5]
chr <- sample(1:19, 4)

# p.mu.chr.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_autosome_counts.csv.gz" # P(mu) 7mer file
# p.mu.chr.methylated.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_autosome_counts.csv.gz"
# p.mu.chr.unmethylated.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_autosome_counts.csv.gz"
# chr <- 19

### IMPORT

k11_chr <- fread(paste0("gunzip -cq ", p.mu.chr.file))
# tmp <- k11_chr
# tmp$chromosome <- 18
# k11_chr$chromosome <- 19
# k11_chr <- rbind(k11_chr, tmp)
# k11_chr <- k11_chr[,c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome")]
# rm(tmp)

k11_chr_meth <- fread(paste0("gunzip -cq ", p.mu.chr.methylated.file))
# tmp <- k11_chr_meth
# tmp$chromosome <- 18
# k11_chr_meth$chromosome <- 19
# k11_chr_meth <- rbind(k11_chr_meth, tmp)
# k11_chr_meth <- k11_chr_meth[,c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome")]
# rm(tmp)

k11_chr_unmeth <- fread(paste0("gunzip -cq ", p.mu.chr.unmethylated.file))
# tmp <- k11_chr_unmeth
# tmp$chromosome <- 18
# k11_chr_unmeth$chromosome <- 19
# k11_chr_unmeth <- rbind(k11_chr_unmeth, tmp)
# k11_chr_unmeth <- k11_chr_unmeth[,c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome")]
# rm(tmp)

### FORMAT

format.fun <- function(k11_chr){
  
colnames(k11_chr) <- c("k11_from", "k11_from_N", "k11_to", "k11_to_N", "chromosome") # colnames
k11_chr$to <- stri_sub(k11_chr$k11_to, 6, -6) # extract middle base changes

## split by chr for training and testing datasets
k11_test <- k11_chr[chromosome %in% chr][,c("k11_from", "k11_from_N", "to", "k11_to_N")] 
k11_train <- k11_chr[!chromosome %in% chr][,c("k11_from", "k11_from_N", "to", "k11_to_N")]
k11_test <- setDT(k11_test)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k11_from, to)] # sum across kmers by chromosome
k11_train <- setDT(k11_train)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k11_from, to)] # sum across kmers by chromosome
rm(k11_chr)

## cut to 7mer 
k7_train <- k11_train # 9mer
colnames(k7_train) <- gsub("k11", "k7", colnames(k11_train)) # colnames
k7_train$k7_from <- stri_sub(k7_train$k7_from, 3, -3)
k7_train <- setDT(k7_train)[, .(k7_from_N = sum(k7_from_N), k7_to_N = sum(k7_to_N)), by = .(k7_from, to)] # sum across kmers

k7_test <- k11_test # 9mer
colnames(k7_test) <- gsub("k11", "k7", colnames(k11_test)) # colnames
k7_test$k7_from <- stri_sub(k7_test$k7_from, 3, -3)
k7_test <- setDT(k7_test)[, .(k7_from_N = sum(k7_from_N), k7_to_N = sum(k7_to_N)), by = .(k7_from, to)] # sum across kmers
rm(k11_test, k11_train)

## average complementary strands
k7_test <- complement(k7_test)
colnames(k7_test) <- c("k7_from", "to", "k7_from_N", "k7_to_N")
k7_train <- complement(k7_train)
colnames(k7_train) <- c("k7_from", "to", "k7_from_N", "k7_to_N")

## subset CG
k7_train <- k7_train[which(stri_sub(k7_train$k7_from, 4, -3) == "CG"),]
k7_test <- k7_test[which(stri_sub(k7_test$k7_from, 4, -3) == "CG"),]
k7_test$k1_from <- stri_sub(k7_test$k7_from, 4, -4) # extract middle base change
k7_train$k1_from <- stri_sub(k7_train$k7_from, 4, -4) # extract middle base change

## calculate pMU for train

k5_train <- k7_train # 5mer
colnames(k5_train) <- gsub("k7", "k5", colnames(k7_train)) # colnames
k5_train$k5_from <- stri_sub(k5_train$k5_from, 2, -2)
k5_train <- setDT(k5_train)[, .(k5_from_N = sum(k5_from_N), k5_to_N = sum(k5_to_N)), by = .(k5_from, to)] # sum across kmers

k3_train <- k7_train # 3mer
colnames(k3_train) <- gsub("k7", "k3", colnames(k7_train))
k3_train$k3_from <- stri_sub(k3_train$k3_from, 3, -3)
k3_train <- setDT(k3_train)[, .(k3_from_N = sum(k3_from_N), k3_to_N = sum(k3_to_N)), by = .(k3_from, to)] # sum across kmers

k1CG_train <- k7_train # k1CG
colnames(k1CG_train) <- gsub("k7", "k1CG", colnames(k7_train))
k1CG_train$k1CG_from <- stri_sub(k1CG_train$k1CG_from, 4, -4)
k1CG_train <- setDT(k1CG_train)[, .(k1CG_from_N = sum(k1CG_from_N), k1CG_to_N = sum(k1CG_to_N)), by = .(k1CG_from, to)] # sum across kmers

# kmer rates of change
k7_train$k7_mu_rate <- k7_train$k7_to_N / k7_train$k7_from_N
k5_train$k5_mu_rate <- k5_train$k5_to_N / k5_train$k5_from_N
k3_train$k3_mu_rate <- k3_train$k3_to_N / k3_train$k3_from_N
k1CG_train$k1CG_mu_rate <- k1CG_train$k1CG_to_N / k1CG_train$k1CG_from_N

k7_train$k7_mu_rate[k7_train$k7_from_N == 0] <- 0
k5_train$k5_mu_rate[k5_train$k5_from_N == 0] <- 0
k3_train$k3_mu_rate[k3_train$k3_from_N == 0] <- 0
k1CG_train$k1CG_mu_rate[k1CG_train$k1CG_from_N == 0] <- 0

k7_train$k5_from <- stri_sub(k7_train$k7_from, 2, -2)
k7_train$k3_from <- stri_sub(k7_train$k7_from, 3, -3)
k7_train$k1_from <- stri_sub(k7_train$k7_from, 4, -4)
k7_train <- k7_train[k7_train, on = c("k7_from", "to")]
k7_train <- k5_train[k7_train, on = c("k5_from", "to")]
k7_train <- k3_train[k7_train, on = c("k3_from", "to")]
k7_train <- k1CG_train[k7_train, on = c("to")]
rm(k3_train, k5_train, k1CG_train)

k7_test$k7_mu_rate <- k7_test$k7_to_N / k7_test$k7_from_N
k7_test$k7_mu_rate[k7_test$k7_from_N == 0] <- 0
k7_test$k1_from <- stri_sub(k7_test$k7_from, 4, -4)

## remove all k7 that are not observed in test
# length(test$k7_from_N[test$k7_from_N == 0]) / nrow(test) # fraction of k7_from_N = 0
# length(train$k7_from_N[train$k7_from_N == 0]) / nrow(train)
ind <- unique(c(which(k7_test$k7_from_N == 0), which(k7_train$k7_from_N == 0)))
if (length(ind) != 0){
  k7_test <- k7_test[-ind,]
  k7_train <- k7_train[-ind,]
} else { 
  k7_test <- k7_test
  k7_train <- k7_train }
# length(test_specific$k7_from_N[test_specific$k7_from_N == 0]) / nrow(test_specific) # fraction of k7_from_N = 0
# length(train_specific$k7_from_N[train_specific$k7_from_N == 0]) / nrow(train_specific)

return(list(k7_train, k7_test))
}

tt_list <- format.fun(k11_chr)
tt_meth_list <- format.fun(k11_chr_meth)
tt_unmeth_list <- format.fun(k11_chr_unmeth)

k7_train <- tt_list[[1]]
k7_test <- tt_list[[2]]

k7_train_meth <- tt_meth_list[[1]]
k7_test_meth <- tt_meth_list[[2]]

k7_train_unmeth <- tt_unmeth_list[[1]]
k7_test_unmeth <- tt_unmeth_list[[2]]

### OBS EXP 

obs.exp <- function(k7_train, k7_test, chr){

out_list <- list()
for (j in 1:length(unique(k7_train$to))){
  out_list[[j]] <- data.table(
    k1_from = k7_train$k1_from[1],
    to = unique(k7_train$to)[j],
    obs_mu = sum(k7_test$k7_mu_rate[k7_test$to == unique(k7_train$to)[j]]) / length(k7_test$k7_mu_rate[k7_test$to == unique(k7_train$to)[j]]), # sum (mu rates) / n (mu rates)
    exp_mu_k1CG = sum(k7_train$k1CG_mu_rate[k7_train$to == unique(k7_train$to)[j]]) / length(k7_train$k1CG_mu_rate[k7_train$to == unique(k7_train$to)[j]]),
    exp_mu_k3 = sum(k7_train$k3_mu_rate[k7_train$to == unique(k7_train$to)[j]]) / length(k7_train$k3_mu_rate[k7_train$to == unique(k7_train$to)[j]]),
    exp_mu_k5 = sum(k7_train$k5_mu_rate[k7_train$to == unique(k7_train$to)[j]]) / length(k7_train$k5_mu_rate[k7_train$to == unique(k7_train$to)[j]]),
    exp_mu_k7 = sum(k7_train$k7_mu_rate[k7_train$to == unique(k7_train$to)[j]]) / length(k7_train$k7_mu_rate[k7_train$to == unique(k7_train$to)[j]])
  )
}
output <- do.call("rbind", out_list)
output$test_chromosome <- paste(chr, collapse = "-")

return(output)
}

meth_A <- obs.exp(k7_train_meth, k7_test_meth, chr)
meth_B <- obs.exp(k7_train, k7_test_meth, chr)
meth_A$meth_status <- "methylated"
meth_B$meth_status <- "none"
out_meth <- rbind(meth_A, meth_B)

unmeth_A <- obs.exp(k7_train_unmeth, k7_test_unmeth, chr)
unmeth_B <- obs.exp(k7_train, k7_test_unmeth, chr)
unmeth_A$meth_status <- "methylated"
unmeth_B$meth_status <- "none"
out_unmeth <- rbind(unmeth_A, unmeth_B)


### EXPORT
fwrite(out_meth, out.specific.meth.file, append = T)
fwrite(out_unmeth, out.specific.unmeth.file, append = T)


#####











