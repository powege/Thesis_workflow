# SCRIPT that creates data table of variables by window for each chromosome

rm(list = ls())
graphics.off()

library(data.table)
# library(evobiR) ##### NOT ON CLUSTER!!!
library(zoo)
library(MASS)
library(stringr)
library(stringi)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {  stop("Arguments must be supplied", call.=FALSE) } 

## set args 
chromo <- as.integer(args[1])
window.size <- as.integer(args[2])
window.shift <- as.integer(args[3])
ref.file <- args[4]
vcf.file <- args[5]
sm.file <- args[6]
Nm.file <- args[7]
cm.file <- args[8]
k7.psnv.file <- args[9]
meth.file <- args[10]
out.file <- args[11]

# chromo <- 19
# window.size <- 1000
# window.shift <- 100
# ref.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr19.txt.gz"
# vcf.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_snps_PASS_chr19.vcf.gz"
# sm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
# Nm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz"
# cm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr19.bed.gz"
# k7.psnv.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV.csv.gz"
# meth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz"
# out.file <- ""

### FUNCTIONS

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

# FUNCTION that sums the middle numbers in a window
n_central <- function(sequence, window.size, window.shift){
  outer.length <- (window.size - window.shift )/2
  central <- sum(sequence[(outer.length + 1):(outer.length + window.shift)])
  return(central)
}

# FUNCTION that vectorises seq
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# FUNCTION that returns complentary kmers (must have columns k7_from and to)
complement <- function(forward){
  
  require(stringi)
  
  complement <- forward # reverse strand
  complement$k7_from <- stri_reverse(complement$k7_from)

  # replace all bases with complement
  complement$k7_from <- gsub("A", "B", complement$k7_from)
  complement$k7_from <- gsub("C", "D", complement$k7_from)
  complement$k7_from <- gsub("T", "A", complement$k7_from)
  complement$k7_from <- gsub("G", "C", complement$k7_from)
  complement$k7_from <- gsub("B", "T", complement$k7_from)
  complement$k7_from <- gsub("D", "G", complement$k7_from)
  complement$k1_from <- gsub("A", "T", complement$k1_from)
  complement$k1_from <- gsub("C", "G", complement$k1_from)

  
  return(complement)
}

### IMPORT

# reference genome
ref <- fread(paste0("gunzip -cq ", ref.file), fill = T, header = F)
# SNVs
vcf <- fread(paste0("gunzip -cq ", vcf.file), fill = T, header = F)
# soft_mask
smask <- fread(paste0("gunzip -cq ", sm.file), fill = T, col.names = c("chromosome", "start", "end"))
# N_mask
Nmask <- fread(paste0("gunzip -cq ", Nm.file), fill = T, col.names = c("chromosome", "start", "end"))
# coverage
cmask <- fread(paste0("gunzip -cq ", cm.file), fill = T, col.names = c("chromosome", "start", "end"))
# k7 psnv
psnv <- fread(paste0("gunzip -cq ", k7.psnv.file), fill = T)
# CG methylation
cgs <- fread(paste0("gunzip -cq ", meth.file), fill = T)

  
### FORMAT ==========

# colnames
colnames(ref) <- c("REF")
ref$POS <- 1:nrow(ref)

# subset chromosome
cmask <- cmask[chromosome == chromo]
smask <- smask[chromosome == chromo]
Nmask <- Nmask[chromosome == chromo]
cgs <- cgs[chromosome == chromo]

# unmethylated CGs
CG_unmeth <- cgs[ES_coverage >= 5 & ES_percent_methylated <= 20][,c("chromosome", "start", "end")]
cgs <- cgs[,c("chromosome", "start", "end")]

# k7mer
b1 <- ref$REF[1:(length(ref$REF)-6)]
b2 <- ref$REF[2:(length(ref$REF)-5)]
b3 <- ref$REF[3:(length(ref$REF)-4)]
b4 <- ref$REF[4:(length(ref$REF)-3)]
b5 <- ref$REF[5:(length(ref$REF)-2)]
b6 <- ref$REF[6:(length(ref$REF)-1)]
b7 <- ref$REF[7:(length(ref$REF))]
k7mer <- paste0(b1, b2, b3, b4, b5, b6, b7)
k7mer[grep("N", k7mer)] <- NA
ref$k7_from <- c(rep(NA, 3), k7mer, rep(NA, 3)) # ensure output equal length to input
rm(b1, b2, b3, b4, b5, b6, b7)

# subset REF by min and max POS 
dt <- ref[POS >= min(which(ref$REF != "N")) & POS <= max(which(ref$REF != "N"))]
rm(ref)
  
# identify POS with SNV
vcf <- unique(vcf$V2)
vcf <- vcf[!is.na(vcf)]
dt$SNV <- 0
dt$SNV[which(dt$POS %in% vcf)] <- 1
rm(vcf)
  
# identify POS with softmask
smask <- unlist(seq2(from = smask$start, to = smask$end))
smask <- unique(smask)
smask <- smask[!is.na(smask)]
dt$soft_mask <- 0
dt$soft_mask[which(dt$POS %in% smask)] <- 1
rm(smask)

# identify POS with N mask
Nmask <- unlist(seq2(from = Nmask$start, to = Nmask$end))
Nmask <- unique(Nmask)
Nmask <- Nmask[!is.na(Nmask)]
dt$N_mask <- 0
dt$N_mask[which(dt$POS %in% Nmask)] <- 1
rm(Nmask)

# identify POS with low coverage
cmask <- unlist(seq2(from = cmask$start, to = cmask$end))
cmask <- unique(cmask)
cmask <- cmask[!is.na(cmask)]
dt$dp_mask <- 0
dt$dp_mask[which(dt$POS %in% cmask)] <- 1
rm(cmask)

# identify POS with CG
cgs <- unlist(seq2(from = cgs$start, to = cgs$end))
cgs <- unique(cgs)
cgs <- cgs[!is.na(cgs)]
dt$CG <- 0
dt$CG[which(dt$POS %in% cgs)] <- 1
rm(cgs)

# identify POS with unmethylated CG
CG_unmeth <- unlist(seq2(from = CG_unmeth$start, to = CG_unmeth$end))
CG_unmeth <- unique(CG_unmeth)
CG_unmeth <- CG_unmeth[!is.na(CG_unmeth)]
dt$CG_unmethylated <- 0
dt$CG_unmethylated[which(dt$POS %in% CG_unmeth)] <- 1
rm(CG_unmeth)

# add psnv
psnv <- psnv[k1_to == "*"] # subset any snv
psnv <- rbind(psnv, complement(psnv)) # expand to include complementary kmers
psnv <- psnv[,c("k7_from", "k7_mu_rate")]
dt <-  psnv[dt, on = c("k7_from")]
rm(psnv)

# calculate variables by sliding window using zoo
x.out <- data.table(
  n_SNV_weighted = rollapply(dt$SNV, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  k7_psnv_weighted = rollapply(dt$k7_mu_rate, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  n_CG_weighted = rollapply(dt$CG, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  n_CG_unmeth_weighted = rollapply(dt$CG_unmethylated, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  n_Nm_weighted = rollapply(dt$N_mask, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  n_sm_weighted = rollapply(dt$soft_mask, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  n_sm_central = rollapply(dt$soft_mask, width = window.size, by = window.shift, FUN = n_central, window.size = window.size, window.shift = window.shift, align = "left"),
  n_cm_weighted = rollapply(dt$dp_mask, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
  n_cm_central = rollapply(dt$dp_mask, width = window.size, by = window.shift, FUN = n_central, window.size = window.size, window.shift = window.shift, align = "left")
)

# Add chromosome, POS from, POS to for each window
outer.length <- (window.size - window.shift )/2
x.out$chromosome <- chromo
x.out$start <- seq(from = (dt$POS[1] + outer.length),
                       to = (dt$POS[1] + (outer.length - 1)) + (window.shift * nrow(x.out)),
                       by = window.shift)
x.out$end <- seq(from <- (dt$POS[1] + (outer.length + window.shift - 1)),
                     to = (dt$POS[1] + (outer.length + window.shift - 1)) + (window.shift * (nrow(x.out)-1)),
                     by = window.shift)

x.out <- x.out[,c("chromosome", "start", "end", 
                  "n_SNV_weighted",
                  "k7_psnv_weighted",
                  "n_CG_weighted",
                  "n_CG_unmeth_weighted",
                  "n_Nm_weighted",
                  "n_sm_weighted",
                  "n_cm_weighted",
                  "n_sm_central",
                  "n_cm_central")]
  
print(paste0(chromo, " done!"))

  
### EXPORT  
fwrite(x.out, out.file, compress = "gzip")


#####
