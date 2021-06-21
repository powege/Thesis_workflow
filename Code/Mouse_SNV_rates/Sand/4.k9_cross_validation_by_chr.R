### SCRIPT that runs cross validation by chromosome (selects 4 random chromosomes for the test dataset)

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

k9.counts.chr.infile <- args[1]
out.file <- args[2]
chr <- as.integer(args[3])

# k9.counts.chr.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k9_SNV_counts_by_chr.csv.gz"
# out.file <- ""
# chr <- 1

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
k9_chr <- fread(paste0("gunzip -cq ", k9.counts.chr.infile))

### FORMAT

## split by chr for training and testing datasets
k9_test <- k9_chr[chromosome %in% chr][,c("k9_from", "k9_to", "k9_from_N", "k9_to_N")] 
k9_train <- k9_chr[!chromosome %in% chr][,c("k9_from", "k9_to", "k9_from_N", "k9_to_N")]
k9_test <- setDT(k9_test)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, k9_to)] # sum across kmers by chromosome
k9_train <- setDT(k9_train)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, k9_to)] # sum across kmers by chromosome
rm(k9_chr)

# sum across complementary kmers
k9_test <- complement(k9_test)
k9_train <- complement(k9_train)

# k1 change
k9_train$k1_from <- stri_sub(k9_train$k9_from, 5, -5)
k9_train$k1_to <- stri_sub(k9_train$k9_to, 5, -5)
k9_test$k1_from <- stri_sub(k9_test$k9_from, 5, -5)
k9_test$k1_to <- stri_sub(k9_test$k9_to, 5, -5)

# subset A and C mutations
k9_test <- k9_test[k1_from %in% c("A", "C")]
k9_train <- k9_train[k1_from %in% c("A", "C")]

## calculate pMU for train

# 9mer
k9_specific <- k9_train

# 7mer
k7_specific <- k9_train
colnames(k7_specific) <- gsub("k9", "k7", colnames(k9_train)) # colnames
k7_specific$k7_from <- stri_sub(k7_specific$k7_from, 2, -2)
k7_specific <- setDT(k7_specific)[, .(k7_from_N = sum(k7_from_N), k7_to_N = sum(k7_to_N)), by = .(k7_from, k1_from, k1_to)] # sum across kmers

# 5mer
k5_specific <- k9_train
colnames(k5_specific) <- gsub("k9", "k5", colnames(k9_train)) # colnames
k5_specific$k5_from <- stri_sub(k5_specific$k5_from, 3, -3)
k5_specific <- setDT(k5_specific)[, .(k5_from_N = sum(k5_from_N), k5_to_N = sum(k5_to_N)), by = .(k5_from, k1_from, k1_to)] # sum across kmers

# 3mer
k3_specific <- k9_train
colnames(k3_specific) <- gsub("k9", "k3", colnames(k9_train))
k3_specific$k3_from <- stri_sub(k3_specific$k3_from, 4, -4)
k3_specific <- setDT(k3_specific)[, .(k3_from_N = sum(k3_from_N), k3_to_N = sum(k3_to_N)), by = .(k3_from, k1_from, k1_to)] # sum across kmers

# 1mer
k1_specific <- k9_train[,c("k1_from", "k1_to", "k9_from_N", "k9_to_N")]
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
k1CG_1 <- k1CG_1[,c("k3_from", "k1_from", "k1_to", "k1CG_group")]
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

k9_train <- as.data.table(rbind.fill(k9_specific, k9_any))
k7_train <- as.data.table(rbind.fill(k7_specific, k7_any))
k5_train <- as.data.table(rbind.fill(k5_specific, k5_any))
k3_train <- as.data.table(rbind.fill(k3_specific, k3_any))
k1_train <- as.data.table(rbind.fill(k1_specific, k1_any))
k1CG_train <- as.data.table(rbind.fill(k1CG_specific, k1CG_any))
rm(k9_specific, k9_any, k7_specific, k7_any, k5_specific, k5_any, k3_specific, k3_any, k1_specific, k1_any, k1CG_specific, k1CG_any)

# kmer rates of change
k9_train$k9_mu_rate <- k9_train$k9_to_N / k9_train$k9_from_N
k7_train$k7_mu_rate <- k7_train$k7_to_N / k7_train$k7_from_N
k5_train$k5_mu_rate <- k5_train$k5_to_N / k5_train$k5_from_N
k3_train$k3_mu_rate <- k3_train$k3_to_N / k3_train$k3_from_N
k1_train$k1_mu_rate <- k1_train$k1_to_N / k1_train$k1_from_N
k1CG_train$k1CG_mu_rate <- k1CG_train$k1CG_to_N / k1CG_train$k1CG_from_N

k9_train$k9_mu_rate[k9_train$k9_from_N == 0] <- 0
k9_train$k9_mu_rate[k9_train$k9_from_N == 0] <- 0
k7_train$k7_mu_rate[k7_train$k7_from_N == 0] <- 0
k5_train$k5_mu_rate[k5_train$k5_from_N == 0] <- 0
k3_train$k3_mu_rate[k3_train$k3_from_N == 0] <- 0
k1_train$k1_mu_rate[k1_train$k1_from_N == 0] <- 0
k1CG_train$k1CG_mu_rate[k1CG_train$k1CG_from_N == 0] <- 0

# mrege
k9_train$k7_from <- stri_sub(k9_train$k9_from, 2, -2)
k9_train$k5_from <- stri_sub(k9_train$k9_from, 3, -3)
k9_train$k3_from <- stri_sub(k9_train$k9_from, 4, -4)
k9_train$k1_from <- stri_sub(k9_train$k9_from, 5, -5)
k9_train <- k7_train[k9_train, on = c("k7_from", "k1_from", "k1_to")]
k9_train <- k5_train[k9_train, on = c("k5_from", "k1_from",  "k1_to")]
k9_train <- k3_train[k9_train, on = c("k3_from", "k1_from", "k1_to")]
k9_train <- k1_train[k9_train, on = c("k1_from", "k1_from", "k1_to")]
k9_train <- k1CG_train[k9_train, on = c("k3_from", "k1_from", "k1_to")]
rm(k1_train, k3_train, k5_train, k7_train, k1CG_train)

## calculate pmu for test
k9_any <- setDT(k9_test)[, .(k9_to_N = sum(k9_to_N)), by = .(k9_from, k1_from, k9_from_N)] # sum across kmers
k9_any$k1_to <- "*"
k9_test <- as.data.table(rbind.fill(k9_test, k9_any))
k9_test$k9_mu_rate <- k9_test$k9_to_N / k9_test$k9_from_N
k9_test$k9_mu_rate[k9_test$k9_from_N == 0] <- 0

## merge test and train
names(k9_test)[names(k9_test) == 'k9_mu_rate'] <- 'k9_mu_rate_obs'
k9_test <- k9_test[,c("k9_from", "k1_to", "k9_mu_rate_obs")]
k9_train <- k9_train[,c("k9_from",
                        "k1_from",
                        "k1_to",
                        "k1_mu_rate",
                        "k1CG_mu_rate",
                        "k3_mu_rate",
                        "k5_mu_rate",
                        "k7_mu_rate",
                        "k9_mu_rate")]
k9_dt <- k9_test[k9_train, on = c("k9_from", "k1_to")]

## calculate absolute error between train and test for each k11
k9_dt$k1_AE <- abs(k9_dt$k1_mu_rate - k9_dt$k9_mu_rate_obs)
k9_dt$k1CG_AE <- abs(k9_dt$k1CG_mu_rate - k9_dt$k9_mu_rate_obs)
k9_dt$k3_AE <- abs(k9_dt$k3_mu_rate - k9_dt$k9_mu_rate_obs)
k9_dt$k5_AE <- abs(k9_dt$k5_mu_rate - k9_dt$k9_mu_rate_obs)
k9_dt$k7_AE <- abs(k9_dt$k7_mu_rate - k9_dt$k9_mu_rate_obs)
k9_dt$k9_AE <- abs(k9_dt$k9_mu_rate - k9_dt$k9_mu_rate_obs)

## calculate the total absolute error for each mutation type 
k9_dt_A <- k9_dt[k1_from == "A"]
out_A <- list()
for (j in 1:length(unique(k9_dt_A$k1_to))){
  out_A[[j]] <- data.table(
    k1_from = k9_dt_A$k1_from[1],
    k1_to = unique(k9_dt_A$k1_to)[j],
    k1_AE = sum(k9_dt_A$k1_AE[k9_dt_A$k1_to == unique(k9_dt_A$k1_to)[j]]),
    k1CG_AE = sum(k9_dt_A$k1CG_AE[k9_dt_A$k1_to == unique(k9_dt_A$k1_to)[j]]),
    k3_AE = sum(k9_dt_A$k3_AE[k9_dt_A$k1_to == unique(k9_dt_A$k1_to)[j]]),
    k5_AE = sum(k9_dt_A$k5_AE[k9_dt_A$k1_to == unique(k9_dt_A$k1_to)[j]]),
    k7_AE = sum(k9_dt_A$k7_AE[k9_dt_A$k1_to == unique(k9_dt_A$k1_to)[j]]),
    k9_AE = sum(k9_dt_A$k9_AE[k9_dt_A$k1_to == unique(k9_dt_A$k1_to)[j]])
    )
}
out_A <- do.call("rbind", out_A)

k9_dt_C <- k9_dt[k1_from == "C"]
out_C <- list()
for (j in 1:length(unique(k9_dt_C$k1_to))){
  out_C[[j]] <- data.table(
    k1_from = k9_dt_C$k1_from[1],
    k1_to = unique(k9_dt_C$k1_to)[j],
    k1_AE = sum(k9_dt_C$k1_AE[k9_dt_C$k1_to == unique(k9_dt_C$k1_to)[j]]),
    k1CG_AE = sum(k9_dt_C$k1CG_AE[k9_dt_C$k1_to == unique(k9_dt_C$k1_to)[j]]),
    k3_AE = sum(k9_dt_C$k3_AE[k9_dt_C$k1_to == unique(k9_dt_C$k1_to)[j]]),
    k5_AE = sum(k9_dt_C$k5_AE[k9_dt_C$k1_to == unique(k9_dt_C$k1_to)[j]]),
    k7_AE = sum(k9_dt_C$k7_AE[k9_dt_C$k1_to == unique(k9_dt_C$k1_to)[j]]),
    k9_AE = sum(k9_dt_C$k9_AE[k9_dt_C$k1_to == unique(k9_dt_C$k1_to)[j]])
    )
}
out_C <- do.call("rbind", out_C)

k9_dt_CG <- k9_dt_C[which(stri_sub(k9_dt_C$k9_from, 5, -4) == "CG"),]
out_CG <- list()
for (j in 1:length(unique(k9_dt_CG$k1_to))){
  out_CG[[j]] <- data.table(
    k1_from = paste0(k9_dt_CG$k1_from[1], "G"),
    k1_to = paste0(unique(k9_dt_CG$k1_to)[j], "G"),
    k1_AE = sum(k9_dt_CG$k1_AE[k9_dt_CG$k1_to == unique(k9_dt_CG$k1_to)[j]]),
    k1CG_AE = sum(k9_dt_CG$k1CG_AE[k9_dt_CG$k1_to == unique(k9_dt_CG$k1_to)[j]]),
    k3_AE = sum(k9_dt_CG$k3_AE[k9_dt_CG$k1_to == unique(k9_dt_CG$k1_to)[j]]),
    k5_AE = sum(k9_dt_CG$k5_AE[k9_dt_CG$k1_to == unique(k9_dt_CG$k1_to)[j]]),
    k7_AE = sum(k9_dt_CG$k7_AE[k9_dt_CG$k1_to == unique(k9_dt_CG$k1_to)[j]]),
    k9_AE = sum(k9_dt_CG$k9_AE[k9_dt_CG$k1_to == unique(k9_dt_CG$k1_to)[j]])
    )
}
out_CG <- do.call("rbind", out_CG)

output <- rbind(out_A, out_C, out_CG)
output$test_chromosome <- paste(chr, collapse = "-")

### EXPORT
fwrite(output, out.file, append = T)


#####
