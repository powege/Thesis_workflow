### Script that calacultes the number of 7mers by annotation 

rm(list=ls())
graphics.off()

library(data.table)
library(stringi)

# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("arguments must be supplied", call.=FALSE)
# } 
# 
# ### Set variables
# 
# ref.file <- args[1]
# sm.file <- args[2]
# cm.file <- args[3]
# gerp.file <- args[4]
# ann.file <- args[5]
# out.file.all.AC <- args[6]
# chr <- as.integer(args[7])

ref.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr19.txt.gz"
sm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
cm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr19.bed.gz"
gerp.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/gerp_constrained_elements.mus_musculus.bed.gz"
ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
meth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz"
chr <- 19

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
kmer_changes <- function(kmer){
  
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
sm <- fread(paste0("gunzip -cq ", sm.file))
gerp <- fread(paste0("gunzip -cq ", gerp.file))
cm <- fread(paste0("gunzip -cq ", cm.file))
ann <- fread(paste0("gunzip -cq ", ann.file))


### FORMAT

# create bed file with annotations of interest (chromosome, start, end, annotation)
# create mask file of soft mask and coverage
# identify kmers
# for each annotation:
  # subset positions and remove masked
  # count 7mers

# annotation file
ann <- ann[chromosome == chr & category %in% c("Exon - CDS",
                                               "Exon - UTR",
                                               "Exon - other",
                                               "Promoter",
                                               "Enhancer - proximal",
                                               "Enhancer - distal",
                                               "CTCF binding",
                                               "Miscellaneous",
                                               "Intron - proximal")][,c("chromosome", "start", "end", "category")]
ann$category[ann$category == "Intron - proximal"] <- "Exon - CDS"
colnames(gerp) <- c("chromosome", "start", "end")
gerp <- gerp[chromosome == chr][, category :=  "GERP"]
all <- rbind(ann, gerp)
all$category <- "All"
annotations <- rbind(ann, gerp, all)
rm(ann, gerp, all)

# mask file
colnames(sm) <- c("chromosome", "start", "end")
colnames(cm) <- c("chromosome", "start", "end")
masked <- rbind(sm, cm)
masked <- masked[chromosome == chr]
rm(sm, cm)

# alignments to continuous reference
colnames(ref) <- c("REF")
ref$POS <- 1:nrow(ref)

# calculate kmers

# k7 ancestral
b1 <- ref$POS[1:(length(ref$POS)-6)]
b2 <- ref$POS[2:(length(ref$POS)-5)]
b3 <- ref$POS[3:(length(ref$POS)-4)]
b4 <- ref$POS[4:(length(ref$POS)-3)]
b5 <- ref$POS[5:(length(ref$POS)-2)]
b6 <- ref$POS[6:(length(ref$POS)-1)]
b7 <- ref$POS[7:(length(ref$POS))]
# k7mer <- paste0(b1, b2, b3, b4, b5, b6, b7)
k7mer <- stringi::stri_c(b1, b2, b3, b4, b5, b6, b7, sep = "")
k7mer[grep("N", k7mer)] <- NA
ref$ANCESTRAL_7MER <- c(rep(NA, 3), k7mer, rep(NA, 3)) # ensure output equal length to input
rm(b1, b2, b3, b4, b5, b6, b7, k7mer)

# cut masked regions
int <- unique(unlist(seq2(from = masked$start, to = masked$end)))
am <- am[!int,]
am <- am[complete.cases(am),]
rm(int, masked)

#####

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
k_out_all$k11_to <- stri_sub(k_out_all$k11_to, 6, -6)
k_out_all_AC <- complement(k_out_all)
k_out_all_AC <- k_out_all_AC[c(which(stri_sub(k_out_all_AC$k11_from, 6, -6) == "A"),
                               which(stri_sub(k_out_all_AC$k11_from, 6, -6) == "C")), ]
k_out_all$chromosome <- chr # add chomosome
k_out_all_AC$chromosome <- chr # add chomosome
rm(k11_all, k11mer_from, mu.count)

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
