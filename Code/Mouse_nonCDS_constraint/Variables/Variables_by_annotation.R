rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS ==========

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

# FUNCTION that vectorises seq
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### SET VARS =========

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {  stop("Arguments must be supplied", call.=FALSE) } 

chromo <- as.integer(args[1])
ann.file <- args[2]
ref.file <- args[3]
vcf.file <- args[4]
sm.file <- args[5]
Nm.file <- args[6]
cm.file <- args[7]
k7.psnv.file <- args[8]
meth.file <- args[9]
out.file <- args[10]

# chromo <- 19
# ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
# ref.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr19.txt.gz"
# vcf.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_snps_PASS_chr19.vcf.gz"
# sm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
# Nm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz"
# cm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr19.bed.gz"
# k7.psnv.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV.csv.gz"
# meth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz"
# out.file <- ""

### IMPORT ==========

# annotations
ann <- fread(paste0("gunzip -cq ", ann.file), fill = T)
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

# subset annotations
ann <- ann[chromosome == chromo & category %in% c("Exon - 5'UTR", 
                                                  "Exon - 3'UTR", 
                                                  "Exon - other",
                                                  "Promoter", 
                                                  "Enhancer - proximal", 
                                                  "Enhancer - distal", 
                                                  "CTCF binding",
                                                  "Miscellaneous", 
                                                  "TAD boundry")][, c("chromosome", 
                                                                      "start", 
                                                                      "end", 
                                                                      "category",  
                                                                      "category_id",  
                                                                      "regulatory_feature_stable_id", 
                                                                      "transcript_id")]
# format reference
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

# add psnv
psnv <- psnv[k1_to == "*"] # subset any snv
psnv <- rbind(psnv, complement(psnv)) # expand to include complementary kmers
psnv <- psnv[,c("k7_from", "k7_mu_rate")]
ref <-  psnv[ref, on = c("k7_from")]
rm(psnv)

# identify POS with SNV
vcf <- unique(vcf$V2)
vcf <- vcf[!is.na(vcf)]
ref$SNV <- 0
ref$SNV[which(ref$POS %in% vcf)] <- 1
rm(vcf)

# identify POS with softmask
smask <- unlist(seq2(from = smask$start, to = smask$end))
smask <- unique(smask)
smask <- smask[!is.na(smask)]
ref$soft_mask <- 0
ref$soft_mask[which(ref$POS %in% smask)] <- 1
rm(smask)

# identify POS with N mask
Nmask <- unlist(seq2(from = Nmask$start, to = Nmask$end))
Nmask <- unique(Nmask)
Nmask <- Nmask[!is.na(Nmask)]
ref$N_mask <- 0
ref$N_mask[which(ref$POS %in% Nmask)] <- 1
rm(Nmask)

# identify POS with low coverage
cmask <- unlist(seq2(from = cmask$start, to = cmask$end))
cmask <- unique(cmask)
cmask <- cmask[!is.na(cmask)]
ref$dp_mask <- 0
ref$dp_mask[which(ref$POS %in% cmask)] <- 1
rm(cmask)

# identify POS with CG
cgs <- unlist(seq2(from = cgs$start, to = cgs$end))
cgs <- unique(cgs)
cgs <- cgs[!is.na(cgs)]
ref$CG <- 0
ref$CG[which(ref$POS %in% cgs)] <- 1
rm(cgs)

# identify POS with unmethylated CG
CG_unmeth <- unlist(seq2(from = CG_unmeth$start, to = CG_unmeth$end))
CG_unmeth <- unique(CG_unmeth)
CG_unmeth <- CG_unmeth[!is.na(CG_unmeth)]
ref$CG_unmethylated <- 0
ref$CG_unmethylated[which(ref$POS %in% CG_unmeth)] <- 1
rm(CG_unmeth)

n_SNV <- rep(NA, nrow(ann))
k7_psnv <- rep(NA, nrow(ann))
n_CG <- rep(NA, nrow(ann))
n_CG_unmethylated <- rep(NA, nrow(ann))
n_Nm <- rep(NA, nrow(ann))
n_sm <- rep(NA, nrow(ann))
n_cm <- rep(NA, nrow(ann))

for (i in 1:nrow(ann)){
    
    n_SNV[i] <- sum(ref$SNV[ann$start[i]:ann$end[i]])
    k7_psnv[i] <- sum(ref$k7_mu_rate[ann$start[i]:ann$end[i]])
    n_CG[i] <- sum(ref$CG[ann$start[i]:ann$end[i]])
    n_CG_unmethylated[i] <- sum(ref$CG_unmethylated[ann$start[i]:ann$end[i]])
    n_Nm[i] <- sum(ref$N_mask[ann$start[i]:ann$end[i]])
    n_sm[i] <- sum(ref$soft_mask[ann$start[i]:ann$end[i]])
    n_cm[i] <- sum(ref$dp_mask[ann$start[i]:ann$end[i]])
    
    print(i)
}

ann$n_SNV <- n_SNV
ann$k7_psnv <- k7_psnv
ann$n_CG <- n_CG
ann$n_CG_unmethylated <- n_CG_unmethylated
ann$n_Nm <- n_Nm
ann$n_sm <- n_sm
ann$n_cm <- n_cm

### EXPORT ==========

fwrite(ann, out.file, append = T, compress = "gzip")










