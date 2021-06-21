### Script that calacultes the number of k9mers and k9mer substitutions by chromosome excluding masked regions

rm(list=ls())
graphics.off()

library(data.table)
library(stringi)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 

### Set variables

ref.file <- args[1]
species1.file <- args[2]
snvs.file <- args[3]
sm.file <- args[4]
cm.file <- args[5]
gerp.file <- args[6]
ann.file <- args[7]
out.file.all <- args[8]
chr <- as.integer(args[9])

# ref.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr19.txt.gz"
# species1.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/mmus_grcm38_v_mpah_v94_alignment_chr19.long.gz"
# snvs.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_snps_PASS_chr19.vcf.gz"
# sm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
# cm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr19.bed.gz"
# gerp.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/gerp_constrained_elements.mus_musculus.bed.gz"
# ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
# chr <- 19

### FUNCTIONS   

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### function that returns all possible kmers and middle base substitutions
k9mer_changes <- function(){

  tmp <- do.call(CJ, replicate(n = 9, expr = c("A", "C", "G", "T"), FALSE))
  from <- do.call(paste, c(tmp, sep=""))
  tmp$V5 <- "A"
  to.A <- do.call(paste, c(tmp, sep=""))
  tmp$V5 <- "C"
  to.C <- do.call(paste, c(tmp, sep=""))
  tmp$V5 <- "G"
  to.G <- do.call(paste, c(tmp, sep=""))
  tmp$V5 <- "T"
  to.T <- do.call(paste, c(tmp, sep=""))
  from <- rep(from, 4)
  to <- c(to.A, to.C, to.G, to.T)
  out <- data.table(from, to)
  out <- out[-which(out$from == out$to),]
  out <- out[order(out$from),]
  colnames(out) <- paste0("k9_", colnames(out))

  return(out)
}

### IMPORT

ref <- fread(paste0("gunzip -cq ", ref.file), header = F)
species1 <- fread(paste0("gunzip -cq ", species1.file))
snvs <- fread(paste0("gunzip -cq ", snvs.file))

sm <- fread(paste0("gunzip -cq ", sm.file))
gerp <- fread(paste0("gunzip -cq ", gerp.file))
cm <- fread(paste0("gunzip -cq ", cm.file))
ann <- fread(paste0("gunzip -cq ", ann.file))

### FORMAT

# mask file
ann <- ann[category %in% c("Exon - CDS", "Exon - UTR", "Exon - other", "Promoter", "Enhancer - proximal", "Enhancer - distal", "CTCF binding", "Miscellaneous", "Intron - proximal")]
ann <- ann[,c("chromosome", "start", "end")]
colnames(sm) <- c("chromosome", "start", "end")
colnames(cm) <- c("chromosome", "start", "end")
colnames(gerp) <- c("chromosome", "start", "end")
masked <- rbind(sm, cm, gerp, ann)
masked <- masked[chromosome == chr]
rm(sm, cm, gerp, ann)

# alignments to continuous reference
colnames(ref) <- c("REF")
colnames(species1) <- c("POS", "REF", "species1")
ref$POS <- 1:nrow(ref)
align <- species1[ref, on = c("POS", "REF")]
rm(ref, species1)

# SNVs to alignment
snvs <- snvs[,c(2,4,5)]
colnames(snvs) <- c("POS", "REF", "ALT")
snvs <- align[snvs, on = c("POS", "REF")]

# infer ancestral and mutant snvs

# ancestral == reference unless alternate == species1
snvs$ancestral <- snvs$REF
snvs$ancestral[which(snvs$ALT == snvs$species1)] <- snvs$ALT[which(snvs$ALT == snvs$species1)]

# mutant == alternate unless alternate == species1
snvs$mutant <- snvs$ALT
snvs$mutant[which(snvs$ALT == snvs$species1)] <- snvs$REF[which(snvs$ALT == snvs$species1)]
snvs <- snvs[,c("POS", "REF", "ancestral", "mutant")]

# align to continuous refernece

align <- align[, c("POS", "REF")]
am <- snvs[align, on = c("POS", "REF")]
rm(snvs, align)
am$ANCESTRAL <- am$REF
am$ANCESTRAL[which(!is.na(am$ancestral))] <- am$ancestral[which(!is.na(am$ancestral))]
am$MUTANT <- am$REF
am$MUTANT[which(!is.na(am$mutant))] <- am$mutant[which(!is.na(am$mutant))]
am <- am[,c("POS", "ANCESTRAL", "MUTANT")]

# calculate kmers

# k9 ancestral
b1 <- am$ANCESTRAL[1:(length(am$ANCESTRAL)-8)]
b2 <- am$ANCESTRAL[2:(length(am$ANCESTRAL)-7)]
b3 <- am$ANCESTRAL[3:(length(am$ANCESTRAL)-6)]
b4 <- am$ANCESTRAL[4:(length(am$ANCESTRAL)-5)]
b5 <- am$ANCESTRAL[5:(length(am$ANCESTRAL)-4)]
b6 <- am$ANCESTRAL[6:(length(am$ANCESTRAL)-3)]
b7 <- am$ANCESTRAL[7:(length(am$ANCESTRAL)-2)]
b8 <- am$ANCESTRAL[8:(length(am$ANCESTRAL)-1)]
b9 <- am$ANCESTRAL[9:(length(am$ANCESTRAL))]
k9mer <- paste0(b1, b2, b3, b4, b5, b6, b7, b8, b9)
k9mer[grep("N", k9mer)] <- NA
am$ANCESTRAL_9MER <- c(rep(NA, 4), k9mer, rep(NA, 4)) # ensure output equal length to input
rm(b1, b2, b3, b4, b5, b6, b7, b8, b9, k9mer)

# k9 mutant
b1 <- am$ANCESTRAL[1:(length(am$ANCESTRAL)-8)]
b2 <- am$ANCESTRAL[2:(length(am$ANCESTRAL)-7)]
b3 <- am$ANCESTRAL[3:(length(am$ANCESTRAL)-6)]
b4 <- am$ANCESTRAL[4:(length(am$ANCESTRAL)-5)]
b5 <- am$MUTANT[5:(length(am$MUTANT)-4)]
b6 <- am$ANCESTRAL[6:(length(am$ANCESTRAL)-3)]
b7 <- am$ANCESTRAL[7:(length(am$ANCESTRAL)-2)]
b8 <- am$ANCESTRAL[8:(length(am$ANCESTRAL)-1)]
b9 <- am$ANCESTRAL[9:(length(am$ANCESTRAL))]
k9mer <- paste0(b1, b2, b3, b4, b5, b6, b7, b8, b9)
k9mer[grep("N", k9mer)] <- NA
am$MUTANT_9MER <- c(rep(NA, 4), k9mer, rep(NA, 4)) # ensure output equal length to input
rm(b1, b2, b3, b4, b5, b6, b7, b8, b9, k9mer)

# cut masked regions
int <- unique(unlist(seq2(from = masked$start, to = masked$end)))
am <- am[!int,]
am <- am[complete.cases(am),]
rm(int, masked)

### kmer counnts excluding methylation

# count ancestral kmers
k9_mu_all <- k9mer_changes() # Identify all possible 9-mer substitutions
k9_all <- data.table(k9_from = unique(k9_mu_all$k9_from))
k9mer_from <- as.data.table(table(am$ANCESTRAL_9MER))
colnames(k9mer_from) <- c(colnames(k9_all), "k9_from_N")
k9mer_from <- k9mer_from[k9_all, on = "k9_from"]
k9mer_from$k9_from_N[is.na(k9mer_from$k9_from_N)] <- 0

# count kmer changes
SNV.POS <- which(am$ANCESTRAL != am$MUTANT) # identify SNV POS
mu.count <- data.table(kmer.from = am$ANCESTRAL_9MER[SNV.POS],
                       kmer.to = am$MUTANT_9MER[SNV.POS])
mu.count <- setDT(mu.count)[,.N ,.(kmer.from, kmer.to)]
colnames(mu.count) <- c("k9_from", "k9_to", "k9_to_N")
mu.count <- mu.count[mu.count$k9_to_N != 0, ]
mu.count <- mu.count[k9_mu_all, on = c("k9_from", "k9_to")]
mu.count$k9_to_N[is.na(mu.count$k9_to_N)] <- 0
k_out_all <- k9mer_from[mu.count, on = "k9_from"]

k_out_all$chromosome <- chr # add chomosome
k_out_all <- k_out_all[,c("chromosome", "k9_from", "k9_to", "k9_from_N", "k9_to_N")]

### OUTPUT
fwrite(x = k_out_all, file = out.file.all, append = T, compress = "gzip")

#######

#####

# # nrow(k_out_all[k9_to_N == 0]) / nrow(k_out_all)
# k_out_all$k9_to <- stri_sub(k_out_all$k9_to, 5, -5)
# # k_out_all_AC <- complement(k_out_all)
# k_out_all_AC <- k_out_all[c(which(stri_sub(k_out_all$k9_from, 5, -5) == "A"),
#                                which(stri_sub(k_out_all$k9_from, 5, -5) == "C")), ]
# # nrow(k_out_all_AC[k9_to_N == 0]) / nrow(k_out_all_AC)
# k3_specific <- k_out_all_AC
# colnames(k3_specific) <- gsub("k9", "k3", colnames(k3_specific)) # colnames
# k3_specific$k3_from <- stri_sub(k3_specific$k3_from, 4, -4)
# k3_specific <- setDT(k3_specific)[, .(k3_from_N = sum(k3_from_N), k3_to_N = sum(k3_to_N)), by = .(k3_from, k3_to)] # sum across kmers
# k3_specific$k3_mu_rate <- k3_specific$k3_to_N / k3_specific$k3_from_N
# hist(k3_specific$k3_mu_rate[which(stri_sub(k3_specific$k3_from, 2, -1) == "CG" & k3_specific$k3_to == "T")])

#####



