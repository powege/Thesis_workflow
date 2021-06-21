# SCRIPT that creates data table of variables by window for each chromosome

rm(list = ls())
graphics.off()

library(data.table)
# library(evobiR) ##### NOT ON CLUSTER!!!
library(zoo)
library(MASS)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {  stop("Arguments must be supplied", call.=FALSE) } 

# set args
species <- args[1]
chromo <- as.integer(args[2])
window.size <- as.integer(args[3])
window.shift <- as.integer(args[4])
ref.file <- args[5]
vcf.file <- args[6]
smask.file <- args[7]
Nmask.file <- args[8]
coverage.file <- args[9]
out.file <- args[10]

# species <- "human"
# chromo <- 21
# window.size <- 750
# window.shift <- 50
# ref.file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr21.csv"
# vcf.file <- "~/Dropbox/PhD/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_controls_chr21.vcf"
# mask.file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_smPOS_Ensembl_GRC38_v94_chr21.csv"
# coverage.file <- "~/Dropbox/PhD/Data/gnomAD/coverage/gnomad_v2.1_fraction90_depth_less10X_chr21.tsv"
# out.file <- ""

### FUNCTIONS

### FUNCTION that returnes the probability of SNV for each k7mer in a sequence. 
## INPUT: vector of bases
## OUTPUT: vector of 7-mer probabilities of SNV
slender_worm <- function(seq, P_SNV){
  
  kmer <- rep(NA, length(seq)-6)
  
  for (i in 1:(length(seq)-6)){
    kmer[i] <- paste0(seq[i], seq[i+1], seq[i+2], seq[i+3], seq[i+4], seq[i+5], seq[i+6])
  }
  
  # identify indices in P_SNV
  ind <- match(kmer, P_SNV$k7_from, nomatch = NA)
  pSNV <- c(rep(NA, 3), P_SNV$p_any_snp_given_k7[ind], rep(NA, 3))
  
  return(pSNV)
}

# FUNCTION that calculates the number of CpG dinucleotides in a sequence (both strands)
n_CpG <- function(seq){
  string <- paste0(seq, collapse = "")
  CpG_count <- str_count(string, "CG") + str_count(string, "GC")
  return(CpG_count)
}

# FUNCTION that calculates the weighted sum of a window
# the middle numbers are weighted 1
# all other numbers are weighted lit a linear function between 0 and 1 (ie 0 outside the window)
sum_weight <- function(sequence, window.size, window.shift){

  outer.length <- (window.size - window.shift )/2
  
  weighted <- sum((sequence[1:outer.length]) * seq(from = 0, to = 1, length.out = outer.length)) +
                   sum((sequence[(outer.length + 1):(outer.length + window.shift)] * 1)) +
                   sum((sequence[(outer.length + window.shift + 1):(window.size)] * seq(from = 1, to = 0, length.out = outer.length))
  )
  return(weighted)
}

# FUNCTION that calculates the weighted number of CpG dinucleotides in a sequence (both strannds)
n_CpG_weight <- function(sequence, window.size, window.shift){
  
  CG.int <- unlist(str_locate_all(string = paste0(sequence, collapse = ""), pattern = c("CG", "GC")))
  binary <- rep(0, length(sequence))
  binary[CG.int] <- 1
  
  outer.length <- (window.size - window.shift )/2
  
  weighted <- sum((binary[1:outer.length]) * seq(from = 0, to = 1, length.out = outer.length)) +
    sum((binary[(outer.length + 1):(outer.length + window.shift)] * 1)) +
    sum((binary[(outer.length + window.shift + 1):(window.size)] * seq(from = 1, to = 0, length.out = outer.length))
  )
  return(weighted)
}

# FUNCTION that sums the middle numbers in a window
n_central <- function(sequence, window.size, window.shift){
  outer.length <- (window.size - window.shift )/2
  central <- sum(sequence[(outer.length + 1):(outer.length + window.shift)])
  return(central)
}

# FUNCTION that vectorises seq
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### IMPORT

# reference genome
ref <- fread(ref.file, fill = T)
# SNVs
vcf <- fread(vcf.file, fill = T)
# soft_mask
smask <- fread(smask.file, fill = T)
# N_mask
Nmask <- fread(Nmask.file, fill = T)
# coverage
coverage <- fread(coverage.file, fill = T)

  
### FORMAT

colnames(ref) <- c("REF")
ref$POS <- 1:nrow(ref)

# subset MAF
if (species == "human"){ vcf <- vcf[vcf$V19 >= 0.001,] }  

# get min and max POS for refenrecne
min.POS <-  min(which(ref$REF != "N"))
max.POS <-   max(which(ref$REF != "N"))
  
# subset REF by min and max POS 
dt <- ref[POS >= min.POS & POS <= max.POS]
rm(ref)
  
# identify POS with SNV
vcf <- unique(vcf$V2)
vcf <- vcf[!is.na(vcf)]
dt$SNV <- 0
dt$SNV[which(dt$POS %in% vcf)] <- 1
rm(vcf)
  
# identify POS with softmask
smask <- unique(smask$V1)
smask <- smask[!is.na(smask)]
dt$soft_mask <- 0
dt$soft_mask[which(dt$POS %in% smask)] <- 1
rm(smask)

# identify POS with N mask
Nmask <- unique(Nmask$V1)
Nmask <- Nmask[!is.na(Nmask)]
dt$N_mask <- 0
dt$N_mask[which(dt$POS %in% Nmask)] <- 1
rm(Nmask)

# identify POS with low coverage
coverage <- unique(coverage$V1)
coverage <- coverage[!is.na(coverage)]
dt$low_coverage <- 0
dt$low_coverage[which(dt$POS %in% coverage)] <- 1
rm(coverage)

# calculate variables by sliding window using zoo
x.out <- data.table(
  # n_SNV = rollapply(dt$SNV, width = window.size, by = window.shift, FUN = sum, align = "left"),
  n_SNV_weighted = rollapply(dt$SNV, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
    
  # n_CG = rollapply(dt$REF, width = window.size, by = window.shift, FUN = n_CpG, align = "left"),
  n_CG_weighted = rollapply(dt$REF, width = window.size, by = window.shift, FUN = n_CpG_weight, window.size = window.size, window.shift = window.shift, align = "left"),
    
  # n_sm = rollapply(dt$soft_mask, width = window.size, by = window.shift, FUN = sum, align = "left"),
  n_sm_weighted = rollapply(dt$soft_mask, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  
  # n_Nm = rollapply(dt$N_mask, width = window.size, by = window.shift, FUN = sum, align = "left"),
  n_Nm_weighted = rollapply(dt$N_mask, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  
  # n_cm = rollapply(dt$low_coverage, width = window.size, by = window.shift, FUN = sum, align = "left"),
  n_cm_weighted = rollapply(dt$low_coverage, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left")
)

# Add chromosome, POS from, POS to for each window
outer.length <- (window.size - window.shift )/2
x.out$chromosome <- rep(chromo, nrow(x.out))
x.out$start <- seq(from = (dt$POS[1] + outer.length),
                       to = (dt$POS[1] + (outer.length - 1)) + (window.shift * nrow(x.out)),
                       by = window.shift)
x.out$end <- seq(from <- (dt$POS[1] + (outer.length + window.shift - 1)),
                     to = (dt$POS[1] + (outer.length + window.shift - 1)) + (window.shift * (nrow(x.out)-1)),
                     by = window.shift)

x.out <- x.out[,c("chromosome", "start", "end", 
                  # "n_SNV", 
                  "n_SNV_weighted",
                  # "n_CG", 
                  "n_CG_weighted",
                  # "n_sm", 
                  "n_sm_weighted",
                  # "n_Nm", 
                  "n_Nm_weighted",
                  # "n_cm", 
                  "n_cm_weighted")]
  
print(paste0(chromo, " done!"))

  
### EXPORT  
fwrite(x.out, out.file)


#####
