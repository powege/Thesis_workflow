### SCRIPT that runs cross validation by chromosome (selects 4 random chromosomes for the test dataset)

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

k11.counts.chr.infile <- args[1]
out.AC.file <- args[2]
out.CCG.file <- args[3]
chr <- sample(1:19, 4)

# k11.counts.chr.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr.csv.gz"
# out.AC.file <- ""
# out.CCG.file <- ""
# chr <- 1

### IMPORT
k11_chr <- fread(paste0("gunzip -cq ", k11.counts.chr.infile))
# k11_chr <- k11_chr[chromosome %in% 1:5]

### FORMAT
names(k11_chr)[names(k11_chr) == 'k11_to'] <- 'to'
k11_chr$k1_from <- stri_sub(k11_chr$k11_from, 6, -6) # extract middle base change

## split by chr for training and testing datasets
k11_test <- k11_chr[chromosome %in% chr][,c("k1_from", "k11_from", "k11_from_N", "to", "k11_to_N")] 
k11_train <- k11_chr[!chromosome %in% chr][,c("k1_from", "k11_from", "k11_from_N", "to", "k11_to_N")]
k11_test <- setDT(k11_test)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k1_from, k11_from, to)] # sum across kmers by chromosome
k11_train <- setDT(k11_train)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k1_from, k11_from, to)] # sum across kmers by chromosome
rm(k11_chr)

## cut to 9mer 
k9_train <- k11_train # 9mer
colnames(k9_train) <- gsub("k11", "k9", colnames(k11_train)) # colnames
k9_train$k9_from <- stri_sub(k9_train$k9_from, 2, -2)
k9_train <- setDT(k9_train)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, to, k1_from)] # sum across kmers

k9_test <- k11_test # 9mer
colnames(k9_test) <- gsub("k11", "k9", colnames(k11_test)) # colnames
k9_test$k9_from <- stri_sub(k9_test$k9_from, 2, -2)
k9_test <- setDT(k9_test)[, .(k9_from_N = sum(k9_from_N), k9_to_N = sum(k9_to_N)), by = .(k9_from, to, k1_from)] # sum across kmers
rm(k11_test, k11_train)

## calculate pMU for train

k7_train <- k9_train # 7mer
colnames(k7_train) <- gsub("k9", "k7", colnames(k9_train)) # colnames
k7_train$k7_from <- stri_sub(k7_train$k7_from, 2, -2)
k7_train <- setDT(k7_train)[, .(k7_from_N = sum(k7_from_N), k7_to_N = sum(k7_to_N)), by = .(k7_from, to)] # sum across kmers

k5_train <- k9_train # 5mer
colnames(k5_train) <- gsub("k9", "k5", colnames(k9_train)) # colnames
k5_train$k5_from <- stri_sub(k5_train$k5_from, 3, -3)
k5_train <- setDT(k5_train)[, .(k5_from_N = sum(k5_from_N), k5_to_N = sum(k5_to_N)), by = .(k5_from, to)] # sum across kmers

k3_train <- k9_train # 3mer
colnames(k3_train) <- gsub("k9", "k3", colnames(k9_train))
k3_train$k3_from <- stri_sub(k3_train$k3_from, 4, -4)
k3_train <- setDT(k3_train)[, .(k3_from_N = sum(k3_from_N), k3_to_N = sum(k3_to_N)), by = .(k3_from, to)] # sum across kmers

k1_train <- k9_train # 3mer
colnames(k1_train) <- gsub("k9", "k1", colnames(k9_train))
k1_train$k1_from <- stri_sub(k1_train$k1_from, 5, -5)
k1_train <- setDT(k1_train)[, .(k1_from_N = sum(k1_from_N), k1_to_N = sum(k1_to_N)), by = .(k1_from, to)] # sum across kmers

k1CG_1 <- k3_train # k1CG
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
k1CG_train <- k1CG_2[k1CG_1, on = c("k1CG_from", "to")]
rm(k1CG_1, k1CG_2)

# kmer rates of change
k9_train$k9_mu_rate <- k9_train$k9_to_N / k9_train$k9_from_N
k7_train$k7_mu_rate <- k7_train$k7_to_N / k7_train$k7_from_N
k5_train$k5_mu_rate <- k5_train$k5_to_N / k5_train$k5_from_N
k3_train$k3_mu_rate <- k3_train$k3_to_N / k3_train$k3_from_N
k1_train$k1_mu_rate <- k1_train$k1_to_N / k1_train$k1_from_N
k1CG_train$k1CG_mu_rate <- k1CG_train$k1CG_to_N / k1CG_train$k1CG_from_N

k9_train$k9_mu_rate[k9_train$k9_from_N == 0] <- 0
k7_train$k7_mu_rate[k7_train$k7_from_N == 0] <- 0
k5_train$k5_mu_rate[k5_train$k5_from_N == 0] <- 0
k3_train$k3_mu_rate[k3_train$k3_from_N == 0] <- 0
k1_train$k1_mu_rate[k1_train$k1_from_N == 0] <- 0
k1CG_train$k1CG_mu_rate[k1CG_train$k1CG_from_N == 0] <- 0

k9_train$k7_from <- stri_sub(k9_train$k9_from, 2, -2)
k9_train$k5_from <- stri_sub(k9_train$k9_from, 3, -3)
k9_train$k3_from <- stri_sub(k9_train$k9_from, 4, -4)
k9_train$k1_from <- stri_sub(k9_train$k9_from, 5, -5)
k9_train <- k7_train[k9_train, on = c("k7_from", "to")]
k9_train <- k5_train[k9_train, on = c("k5_from", "to")]
k9_train <- k3_train[k9_train, on = c("k3_from", "to")]
k9_train <- k1_train[k9_train, on = c("k1_from", "to")]
k9_train <- k1CG_train[k9_train, on = c("k3_from", "to")]
rm(k1_train, k3_train, k5_train, k7_train, k1CG_train)

## calculate pmu for test
k9_test$k9_mu_rate <- k9_test$k9_to_N / k9_test$k9_from_N
k9_test$k9_mu_rate[k9_test$k9_from_N == 0] <- 0
k9_test$k1_from <- stri_sub(k9_test$k9_from, 5, -5)

## merge test and train
names(k9_test)[names(k9_test) == 'k9_mu_rate'] <- 'k9_mu_rate_obs'
k9_test <- k9_test[,c("k9_from", "to", "k9_mu_rate_obs")]
k9_train <- k9_train[,c("k9_from",
                        "k1_from",
                        "to",
                        "k1CG_mu_rate",
                        "k1_mu_rate",
                        "k3_mu_rate",
                        "k5_mu_rate",
                        "k7_mu_rate",
                        "k9_mu_rate")]
k9_dt <- k9_test[k9_train, on = c("k9_from", "to")]

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
for (j in 1:length(unique(k9_dt_A$to))){
  out_A[[j]] <- data.table(
    k1_from = k9_dt_A$k1_from[1],
    to = unique(k9_dt_A$to)[j],
    k1_AE = sum(k9_dt_A$k1_AE[k9_dt_A$to == unique(k9_dt_A$to)[j]]),
    k1CG_AE = sum(k9_dt_A$k1CG_AE[k9_dt_A$to == unique(k9_dt_A$to)[j]]),
    k3_AE = sum(k9_dt_A$k3_AE[k9_dt_A$to == unique(k9_dt_A$to)[j]]),
    k5_AE = sum(k9_dt_A$k5_AE[k9_dt_A$to == unique(k9_dt_A$to)[j]]),
    k7_AE = sum(k9_dt_A$k7_AE[k9_dt_A$to == unique(k9_dt_A$to)[j]]),
    k9_AE = sum(k9_dt_A$k9_AE[k9_dt_A$to == unique(k9_dt_A$to)[j]])
    )
}

k9_dt_C <- k9_dt[k1_from == "C"]
out_C <- list()
for (j in 1:length(unique(k9_dt_C$to))){
  out_C[[j]] <- data.table(
    k1_from = k9_dt_C$k1_from[1],
    to = unique(k9_dt_C$to)[j],
    k1_AE = sum(k9_dt_C$k1_AE[k9_dt_C$to == unique(k9_dt_C$to)[j]]),
    k1CG_AE = sum(k9_dt_C$k1CG_AE[k9_dt_C$to == unique(k9_dt_C$to)[j]]),
    k3_AE = sum(k9_dt_C$k3_AE[k9_dt_C$to == unique(k9_dt_C$to)[j]]),
    k5_AE = sum(k9_dt_C$k5_AE[k9_dt_C$to == unique(k9_dt_C$to)[j]]),
    k7_AE = sum(k9_dt_C$k7_AE[k9_dt_C$to == unique(k9_dt_C$to)[j]]),
    k9_AE = sum(k9_dt_C$k9_AE[k9_dt_C$to == unique(k9_dt_C$to)[j]])
    )
}
out_A_C <- do.call("rbind", c(out_A, out_C))
out_A_C$test_chromosome <- paste(chr, collapse = "-")

k9_dt_C_noCG <- k9_dt_C[which(stri_sub(k9_dt_C$k9_from, 5, -4) != "CG"),]
out_C_noCG <- list()
for (j in 1:length(unique(k9_dt_C_noCG$to))){
  out_C_noCG[[j]] <- data.table(
    k1_from = k9_dt_C_noCG$k1_from[1],
    to = unique(k9_dt_C_noCG$to)[j],
    k1_AE = sum(k9_dt_C_noCG$k1_AE[k9_dt_C_noCG$to == unique(k9_dt_C_noCG$to)[j]]),
    k1CG_AE = sum(k9_dt_C_noCG$k1CG_AE[k9_dt_C_noCG$to == unique(k9_dt_C_noCG$to)[j]]),
    k3_AE = sum(k9_dt_C_noCG$k3_AE[k9_dt_C_noCG$to == unique(k9_dt_C_noCG$to)[j]]),
    k5_AE = sum(k9_dt_C_noCG$k5_AE[k9_dt_C_noCG$to == unique(k9_dt_C_noCG$to)[j]]),
    k7_AE = sum(k9_dt_C_noCG$k7_AE[k9_dt_C_noCG$to == unique(k9_dt_C_noCG$to)[j]]),
    k9_AE = sum(k9_dt_C_noCG$k9_AE[k9_dt_C_noCG$to == unique(k9_dt_C_noCG$to)[j]])
    )
}
out_C_noCG <- do.call("rbind", out_C_noCG)
out_C_noCG$test_chromosome <- paste(chr, collapse = "-")
out_C_noCG$CG_status <- "nonCG"

k9_dt_C_CG <- k9_dt_C[which(stri_sub(k9_dt_C$k9_from, 5, -4) == "CG"),]
out_C_CG <- list()
for (j in 1:length(unique(k9_dt_C_CG$to))){
  out_C_CG[[j]] <- data.table(
    k1_from = k9_dt_C_CG$k1_from[1],
    to = unique(k9_dt_C_CG$to)[j],
    k1_AE = sum(k9_dt_C_CG$k1_AE[k9_dt_C_CG$to == unique(k9_dt_C_CG$to)[j]]),
    k1CG_AE = sum(k9_dt_C_CG$k1CG_AE[k9_dt_C_CG$to == unique(k9_dt_C_CG$to)[j]]),
    k3_AE = sum(k9_dt_C_CG$k3_AE[k9_dt_C_CG$to == unique(k9_dt_C_CG$to)[j]]),
    k5_AE = sum(k9_dt_C_CG$k5_AE[k9_dt_C_CG$to == unique(k9_dt_C_CG$to)[j]]),
    k7_AE = sum(k9_dt_C_CG$k7_AE[k9_dt_C_CG$to == unique(k9_dt_C_CG$to)[j]]),
    k9_AE = sum(k9_dt_C_CG$k9_AE[k9_dt_C_CG$to == unique(k9_dt_C_CG$to)[j]])
    )
}
out_C_CG <- do.call("rbind", out_C_CG)
out_C_CG$test_chromosome <- paste(chr, collapse = "-")
out_C_CG$CG_status <- "CG"
out_C_CG <- rbind(out_C_CG, out_C_noCG)

### EXPORT
fwrite(out_A_C, out.AC.file, append = T)
fwrite(out_C_CG, out.CCG.file, append = T)


#####
