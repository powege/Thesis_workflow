### Script that calacultes the number of k11mers and k11mer substitutions by chromosome excluding masked regions

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
species2.file <- args[3]
snvs.file <- args[4]
sm.file <- args[5]
cm.file <- args[6]
gerp.file <- args[7]
ann.file <- args[8]
meth.file <- args[9]
out.file.all <- args[10]
out.file.CGmeth <- args[11]
out.file.CGunmeth <- args[12]
out.file.all.AC <- args[13]
out.file.CGmeth.C <- args[14]
out.file.CGunmeth.C <- args[15]
chr <- as.integer(args[16])

# ref.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr19.txt.gz"
# species1.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/mmus_grcm38_v_mpah_v94_alignment_chr19.long.gz"
# species2.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/mmus_grcm38_v_mpah_v94_alignment_chr19.long.gz"
# snvs.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_snps_PASS_chr19.vcf.gz"
# sm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
# cm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr19.bed.gz"
# gerp.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/gerp_constrained_elements.mus_musculus.bed.gz"
# ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
# meth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz"
# chr <- 19

# ref.file <- "/well/lindgren/George//Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr19.txt.gz"
# species1.file <- "/well/lindgren/George//Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/mmus_grcm38_v_mpah_v94_alignment_chr19.long.gz"
# species2.file <- "/well/lindgren/George//Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/mmus_grcm38_v_mpah_v94_alignment_chr19.long.gz"
# snvs.file <- "/well/lindgren/George//Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_snps_PASS_chr19.vcf.gz"
# sm.file <- "/well/lindgren/George//Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
# cm.file <- "/well/lindgren/George//Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr19.bed.gz"
# gerp.file <- "/well/lindgren/George//Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/ensembl_v94_gerp_constrained_elements.mus_musculus.bed.gz"
# ann.file <- "/well/lindgren/George//Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
# out.file <- ""
# chr <- 19

### FUNCTIONS   

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### function that returns all possible kmers and middle base substitutions
k11mer_changes <- function(kmer){

  tmp <- do.call(CJ, replicate(n = kmer, expr = c("A", "C", "G", "T"), FALSE))
  from <- do.call(paste, c(tmp, sep=""))
  tmp$V6 <- "A"
  to.A <- do.call(paste, c(tmp, sep=""))
  tmp$V6 <- "C"
  to.C <- do.call(paste, c(tmp, sep=""))
  tmp$V6 <- "G"
  to.G <- do.call(paste, c(tmp, sep=""))
  tmp$V6 <- "T"
  to.T <- do.call(paste, c(tmp, sep=""))
  from <- rep(from, 4)
  to <- c(to.A, to.C, to.G, to.T)
  out <- data.table(from, to)
  out <- out[-which(out$from == out$to),]
  out <- out[order(out$from),]
  colnames(out) <- paste0("k", kmer, "_", colnames(out))

  return(out)
}

# FUNCTION that sums across complimantary strands (ie AGC == GCT)
complement <- function(forward){
  
  complement <- forward # reverse strand
  complement$k11_from <- stri_reverse(complement$k11_from)
  complement$k11_to <- stri_reverse(complement$k11_to)
  
  # replace all bases with complement
  complement$k11_from <- gsub("A", "B", complement$k11_from)
  complement$k11_from <- gsub("C", "D", complement$k11_from)
  complement$k11_from <- gsub("T", "A", complement$k11_from)
  complement$k11_from <- gsub("G", "C", complement$k11_from)
  complement$k11_from <- gsub("B", "T", complement$k11_from)
  complement$k11_from <- gsub("D", "G", complement$k11_from)
  complement$k11_to <- gsub("A", "B", complement$k11_to)
  complement$k11_to <- gsub("C", "D", complement$k11_to)
  complement$k11_to <- gsub("T", "A", complement$k11_to)
  complement$k11_to <- gsub("G", "C", complement$k11_to)
  complement$k11_to <- gsub("B", "T", complement$k11_to)
  complement$k11_to <- gsub("D", "G", complement$k11_to)
  
  # sum across complements
  output <- rbind(forward, complement)
  output <- setDT(output)[, .(k11_from_N = sum(k11_from_N), k11_to_N = sum(k11_to_N)), by = .(k11_from, k11_to)]
  
  return(output)
}

### IMPORT

ref <- fread(paste0("gunzip -cq ", ref.file), header = F)
species1 <- fread(paste0("gunzip -cq ", species1.file))
species2 <- fread(paste0("gunzip -cq ", species2.file))
snvs <- fread(paste0("gunzip -cq ", snvs.file))

sm <- fread(paste0("gunzip -cq ", sm.file))
gerp <- fread(paste0("gunzip -cq ", gerp.file))
cm <- fread(paste0("gunzip -cq ", cm.file))
ann <- fread(paste0("gunzip -cq ", ann.file))

meth_all <- fread(paste0("gunzip -cq ", meth.file))

### FORMAT

# methylation
unmeth <- meth_all[chromosome == chr & coverage > 5 & percent_methylated < 20][,c("chromosome", "start", "end")]
meth <- meth_all[chromosome == chr & coverage > 5 & percent_methylated > 60][,c("chromosome", "start", "end")]
rm(meth_all)

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
colnames(species2) <- c("POS", "REF", "species2")

ref$POS <- 1:nrow(ref)

align <- species1[ref, on = c("POS", "REF")]
align <- species2[align, on = c("POS", "REF")]

rm(ref, species1, species2)

# SNVs to alignment
snvs <- snvs[,c(2,4,5)]
colnames(snvs) <- c("POS", "REF", "ALT")
snvs <- align[snvs, on = c("POS", "REF")]

# infer ancestral and mutant snvs

# ancestral == reference unless alternate == species1, or if no species1 alignment, alternate == species2 
snvs$ancestral <- snvs$REF
snvs$ancestral[which(snvs$ALT == snvs$species1)] <- snvs$ALT[which(snvs$ALT == snvs$species1)]
snvs$ancestral[which(is.na(snvs$species1) & snvs$ALT == snvs$species2)] <- snvs$ALT[which(is.na(snvs$species1) & snvs$ALT == snvs$species2)] 

# mutant == alternate unless alternate == species1, or if no species1 alignment, alternate == species2 
snvs$mutant <- snvs$ALT
snvs$mutant[which(snvs$ALT == snvs$species1)] <- snvs$REF[which(snvs$ALT == snvs$species1)]
snvs$mutant[which(is.na(snvs$species1) & snvs$ALT == snvs$species2)] <- snvs$REF[which(is.na(snvs$species1) & snvs$ALT == snvs$species2)] 
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

# k11 ancestral
b1 <- am$ANCESTRAL[1:(length(am$ANCESTRAL)-10)]
b2 <- am$ANCESTRAL[2:(length(am$ANCESTRAL)-9)]
b3 <- am$ANCESTRAL[3:(length(am$ANCESTRAL)-8)]
b4 <- am$ANCESTRAL[4:(length(am$ANCESTRAL)-7)]
b5 <- am$ANCESTRAL[5:(length(am$ANCESTRAL)-6)]
b6 <- am$ANCESTRAL[6:(length(am$ANCESTRAL)-5)]
b7 <- am$ANCESTRAL[7:(length(am$ANCESTRAL)-4)]
b8 <- am$ANCESTRAL[8:(length(am$ANCESTRAL)-3)]
b9 <- am$ANCESTRAL[9:(length(am$ANCESTRAL)-2)]
b10 <- am$ANCESTRAL[10:(length(am$ANCESTRAL)-1)]
b11 <- am$ANCESTRAL[11:(length(am$ANCESTRAL))]
k11mer <- paste0(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11)
k11mer[grep("N", k11mer)] <- NA
am$ANCESTRAL_11MER <- c(rep(NA, 5), k11mer, rep(NA, 5)) # ensure output equal length to input
rm(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, k11mer)

# k11 mutant
b1 <- am$ANCESTRAL[1:(length(am$ANCESTRAL)-10)]
b2 <- am$ANCESTRAL[2:(length(am$ANCESTRAL)-9)]
b3 <- am$ANCESTRAL[3:(length(am$ANCESTRAL)-8)]
b4 <- am$ANCESTRAL[4:(length(am$ANCESTRAL)-7)]
b5 <- am$ANCESTRAL[5:(length(am$ANCESTRAL)-6)]
b6 <- am$MUTANT[6:(length(am$MUTANT)-5)]
b7 <- am$ANCESTRAL[7:(length(am$ANCESTRAL)-4)]
b8 <- am$ANCESTRAL[8:(length(am$ANCESTRAL)-3)]
b9 <- am$ANCESTRAL[9:(length(am$ANCESTRAL)-2)]
b10 <- am$ANCESTRAL[10:(length(am$ANCESTRAL)-1)]
b11 <- am$ANCESTRAL[11:(length(am$ANCESTRAL))]
k11mer <- paste0(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11)
k11mer[grep("N", k11mer)] <- NA
am$MUTANT_11MER <- c(rep(NA, 5), k11mer, rep(NA, 5)) # ensure output equal length to input
rm(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, k11mer)

# cut masked regions
int <- unique(unlist(seq2(from = masked$start, to = masked$end)))
am <- am[!int,]
am <- am[complete.cases(am),]
rm(int, masked)

### kmer counnts excluding methylation

# count ancestral kmers
k11_mu_all <- k11mer_changes(11) # Identify all possible 11-mer substitutions
k11_all <- data.table(k11_from = unique(k11_mu_all$k11_from))
k11mer_from <- as.data.table(table(am$ANCESTRAL_11MER))
colnames(k11mer_from) <- c(colnames(k11_all), "k11_from_N")
k11mer_from <- k11mer_from[k11_all, on = "k11_from"]
k11mer_from$k11_from_N[is.na(k11mer_from$k11_from_N)] <- 0

# count kmer changes
SNV.POS <- which(am$ANCESTRAL != am$MUTANT) # identify SNV POS
mu.count <- data.table(kmer.from = am$ANCESTRAL_11MER[SNV.POS],
                       kmer.to = am$MUTANT_11MER[SNV.POS])
mu.count <- setDT(mu.count)[,.N ,.(kmer.from, kmer.to)]
colnames(mu.count) <- c("k11_from", "k11_to", "k11_to_N")
mu.count <- mu.count[mu.count$k11_to_N != 0, ]
mu.count <- mu.count[k11_mu_all, on = c("k11_from", "k11_to")]
mu.count$k11_to_N[is.na(mu.count$k11_to_N)] <- 0

k_out_all <- k11mer_from[mu.count, on = "k11_from"]
# nrow(k_out_all[k11_to_N == 0]) / nrow(k_out_all)

k_out_all$k11_to <- stri_sub(k_out_all$k11_to, 6, -6)
k_out_all_AC <- complement(k_out_all)
k_out_all_AC <- k_out_all_AC[c(which(stri_sub(k_out_all_AC$k11_from, 6, -6) == "A"),
                               which(stri_sub(k_out_all_AC$k11_from, 6, -6) == "C")), ]
k_out_all$chromosome <- chr # add chomosome
k_out_all_AC$chromosome <- chr # add chomosome
rm(k11_all, k11mer_from, mu.count)

nrow(k_out_all_AC[k11_to_N == 0]) / nrow(k_out_all_AC)



### kmer counts methylated CG

# count ancestral kmers
k11_mu_all_CG <- k11_mu_all[c(which(stri_sub(k11_mu_all$k11_from, 5, -6) == "CG"),
                    which(stri_sub(k11_mu_all$k11_from, 6, -5) == "CG")), ]
k11_all_CG <- data.table(k11_from = unique(k11_mu_all_CG$k11_from))
am_CG <- am[c(which(stri_sub(am$ANCESTRAL_11MER, 5, -6) == "CG"),
                      which(stri_sub(am$ANCESTRAL_11MER, 6, -5) == "CG")), ]
am_methylated <- am_CG[am_CG$POS %in% c(meth$start, meth$end),]
k11mer_from <- as.data.table(table(am_methylated$ANCESTRAL_11MER))
colnames(k11mer_from) <- c(colnames(k11_all_CG), "k11_from_N")
k11mer_from <- k11mer_from[k11_all_CG, on = "k11_from"]
k11mer_from$k11_from_N[is.na(k11mer_from$k11_from_N)] <- 0

# count kmer changes
SNV.POS <- which(am_methylated$ANCESTRAL != am_methylated$MUTANT) # identify SNV POS
mu.count <- data.table(kmer.from = am_methylated$ANCESTRAL_11MER[SNV.POS],
                       kmer.to = am_methylated$MUTANT_11MER[SNV.POS])
mu.count <- setDT(mu.count)[,.N ,.(kmer.from, kmer.to)]
colnames(mu.count) <- c("k11_from", "k11_to", "k11_to_N")
mu.count <- mu.count[mu.count$k11_to_N != 0, ]
mu.count <- mu.count[k11_mu_all_CG, on = c("k11_from", "k11_to")]
mu.count$k11_to_N[is.na(mu.count$k11_to_N)] <- 0

k_out_meth <- k11mer_from[mu.count, on = "k11_from"]
k_out_meth$k11_to <- stri_sub(k_out_meth$k11_to, 6, -6)
k_out_meth_C <- complement(k_out_meth)
k_out_meth_C <- k_out_meth_C[which(stri_sub(k_out_meth_C$k11_from, 6, -6) == "C"), ]
k_out_meth$chromosome <- chr # add chomosome
k_out_meth_C$chromosome <- chr # add chomosome
rm(k11_mu_all, mu.count)

### kmer counts unmethylated CG
am_unmethylated <- am_CG[am_CG$POS %in% c(unmeth$start, unmeth$end),]
k11mer_from <- as.data.table(table(am_unmethylated$ANCESTRAL_11MER))
colnames(k11mer_from) <- c(colnames(k11_all_CG), "k11_from_N")
k11mer_from <- k11mer_from[k11_all_CG, on = "k11_from"]
k11mer_from$k11_from_N[is.na(k11mer_from$k11_from_N)] <- 0

SNV.POS <- which(am_unmethylated$ANCESTRAL != am_unmethylated$MUTANT) # identify SNV POS
mu.count <- data.table(kmer.from = am_unmethylated$ANCESTRAL_11MER[SNV.POS],
                       kmer.to = am_unmethylated$MUTANT_11MER[SNV.POS])
mu.count <- setDT(mu.count)[,.N ,.(kmer.from, kmer.to)]
colnames(mu.count) <- c("k11_from", "k11_to", "k11_to_N")
mu.count <- mu.count[mu.count$k11_to_N != 0, ]
mu.count <- mu.count[k11_mu_all_CG, on = c("k11_from", "k11_to")]
mu.count$k11_to_N[is.na(mu.count$k11_to_N)] <- 0

k_out_unmeth <- k11mer_from[mu.count, on = "k11_from"]
k_out_unmeth$k11_to <- stri_sub(k_out_unmeth$k11_to, 6, -6)
k_out_unmeth_C <- complement(k_out_unmeth)
k_out_unmeth_C <- k_out_unmeth_C[which(stri_sub(k_out_unmeth_C$k11_from, 6, -6) == "C"), ]
k_out_unmeth_C$chromosome <- chr # add chomosome
k_out_unmeth$chromosome <- chr # add chomosome
rm(k11mer_from, mu.count)


### OUTPUT
fwrite(x = k_out_all, file = out.file.all, append = T, compress = "gzip")
fwrite(x = k_out_meth, file = out.file.CGmeth, append = T, compress = "gzip")
fwrite(x = k_out_unmeth, file = out.file.CGunmeth, append = T, compress = "gzip")
fwrite(x = k_out_all_AC, file = out.file.all.AC, append = T, compress = "gzip")
fwrite(x = k_out_meth_C, file = out.file.CGmeth.C, append = T, compress = "gzip")
fwrite(x = k_out_unmeth_C, file = out.file.CGunmeth.C, append = T, compress = "gzip")

#######
