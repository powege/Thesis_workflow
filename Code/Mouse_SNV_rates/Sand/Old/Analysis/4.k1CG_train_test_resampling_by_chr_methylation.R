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
k11.counts.chr.meth.infile <- args[2]
out.file <- args[3]
chr <- sample(1:19, 4)

# k11.counts.chr.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr.csv.gz"
# k11.counts.chr.meth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr.csv.gz"
# chr <- 1

### IMPORT
k11_chr <- fread(paste0("gunzip -cq ", k11.counts.chr.infile))
# k11_chr <- k11_chr[chromosome %in% 1:5]
k11_chr_meth <- fread(paste0("gunzip -cq ", k11.counts.chr.meth.infile))

### FORMAT
names(k11_chr)[names(k11_chr) == 'k11_to'] <- 'to'
k11_chr$k1_from <- stri_sub(k11_chr$k11_from, 6, -6) # extract middle base change

names(k11_chr_meth)[names(k11_chr_meth) == 'k11_to'] <- 'to'
k11_chr_meth$k1_from <- stri_sub(k11_chr_meth$k11_from, 6, -6) # extract middle base change

## split by chr for training and testing datasets
k11_train <- k11_chr[!chromosome %in% chr][,c("k1_from", "k11_from", "k11_from_N", "to", "k11_to_N")]
k11_train <- setDT(k11_train)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k1_from, k11_from, to)] # sum across kmers by chromosome
rm(k11_chr)

k11_test_meth <- k11_chr_meth[chromosome %in% chr][,c("k1_from", "k11_from", "k11_from_N", "to", "k11_to_N")] 
k11_train_meth <- k11_chr_meth[!chromosome %in% chr][,c("k1_from", "k11_from", "k11_from_N", "to", "k11_to_N")]
k11_test_meth <- setDT(k11_test_meth)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k1_from, k11_from, to)] # sum across kmers by chromosome
k11_train_meth <- setDT(k11_train_meth)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k1_from, k11_from, to)] # sum across kmers by chromosome
rm(k11_chr_meth)

## cut to k1CG
k1CG_train <- k11_train 
k1CG_train <- k1CG_train[which(stri_sub(k1CG_train$k11_from, 6, -5) == "CG"),]
colnames(k1CG_train) <- gsub("k11", "k1CG", colnames(k11_train)) # colnames
k1CG_train$k1CG_from <- stri_sub(k1CG_train$k1CG_from, 6, -6)
k1CG_train <- setDT(k1CG_train)[, .(k1CG_from_N = sum(k1CG_from_N), k1CG_to_N = sum(k1CG_to_N)), by = .(k1CG_from, to, k1_from)] # sum across kmers

k1CG_train_meth <- k11_train_meth 
colnames(k1CG_train_meth) <- gsub("k11", "k1CG", colnames(k11_train_meth)) # colnames
k1CG_train_meth$k1CG_from <- stri_sub(k1CG_train_meth$k1CG_from, 6, -6)
k1CG_train_meth <- setDT(k1CG_train_meth)[, .(k1CG_from_N = sum(k1CG_from_N), k1CG_to_N = sum(k1CG_to_N)), by = .(k1CG_from, to, k1_from)] # sum across kmers

k1CG_test_meth <- k11_test_meth 
colnames(k1CG_test_meth) <- gsub("k11", "k1CG", colnames(k11_test_meth)) # colnames
k1CG_test_meth$k1CG_from <- stri_sub(k1CG_test_meth$k1CG_from, 6, -6)
k1CG_test_meth <- setDT(k1CG_test_meth)[, .(k1CG_from_N = sum(k1CG_from_N), k1CG_to_N = sum(k1CG_to_N)), by = .(k1CG_from, to, k1_from)] # sum across kmers

rm(k11_test_meth, k11_train, k11_train_meth)

# kmer rates of change
k1CG_train$k1CG_mu_rate <- k1CG_train$k1CG_to_N / k1CG_train$k1CG_from_N
k1CG_train_meth$k1CG_mu_rate <- k1CG_train_meth$k1CG_to_N / k1CG_train_meth$k1CG_from_N
k1CG_test_meth$k1CG_mu_rate <- k1CG_test_meth$k1CG_to_N / k1CG_test_meth$k1CG_from_N

k1CG_train$k1CG_mu_rate[k1CG_train$k1CG_from_N == 0] <- 0
k1CG_train_meth$k1CG_mu_rate[k1CG_train_meth$k1CG_from_N == 0] <- 0
k1CG_test_meth$k1CG_mu_rate[k1CG_test_meth$k1CG_from_N == 0] <- 0

## merge test and train
names(k1CG_test_meth)[names(k1CG_test_meth) == 'k1CG_mu_rate'] <- 'k1CG_mu_rate_obs'
names(k1CG_train_meth)[names(k1CG_train_meth) == 'k1CG_mu_rate'] <- 'k1CG_mu_rate_meth'
k1CG_test_meth <- k1CG_test_meth[,c("k1CG_from", "to", "k1CG_mu_rate_obs")]
k1CG_train <- k1CG_train[,c("k1CG_from","to","k1CG_mu_rate")]
k1CG_train_meth <- k1CG_train_meth[,c("k1CG_from","to","k1CG_mu_rate_meth")]
k1CG_dt <- k1CG_test_meth[k1CG_train_meth, on = c("k1CG_from", "to")]
k1CG_dt <- k1CG_dt[k1CG_train, on = c("k1CG_from", "to")]

## calculate absolute error between train and test for each k11
k1CG_dt$k1CG_meth_AE <- abs(k1CG_dt$k1CG_mu_rate_meth - k1CG_dt$k1CG_mu_rate_obs)
k1CG_dt$k1CG_AE <- abs(k1CG_dt$k1CG_mu_rate - k1CG_dt$k1CG_mu_rate_obs)
k1CG_dt$test_chromosome <- paste(chr, collapse = "-")


### EXPORT
fwrite(k1CG_dt, out.file, append = T)


#####
