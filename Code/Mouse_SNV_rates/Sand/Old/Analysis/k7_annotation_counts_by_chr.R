### Script that calacultes the number of k11mers and k11mer substitutions by chromosome
# excluding masked regions

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
sm.file <- args[2]
cm.file <- args[3]
gerp.file <- args[4]
ann.file <- args[5]
out.file <- args[6]
chr <- as.integer(args[7])

# ref.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr19.txt.gz"
# sm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
# cm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr19.bed.gz"
# gerp.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/gerp_constrained_elements.mus_musculus.bed.gz"
# ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
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
  complement$k7_from <- stri_reverse(complement$k7_from)

  # replace all bases with complement
  complement$k7_from <- gsub("A", "B", complement$k7_from)
  complement$k7_from <- gsub("C", "D", complement$k7_from)
  complement$k7_from <- gsub("T", "A", complement$k7_from)
  complement$k7_from <- gsub("G", "C", complement$k7_from)
  complement$k7_from <- gsub("B", "T", complement$k7_from)
  complement$k7_from <- gsub("D", "G", complement$k7_from)
  
  # sum across complements
  output <- rbind(forward, complement)
  output <- setDT(output)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]
  
  return(output)
}

### IMPORT

ref <- fread(paste0("gunzip -cq ", ref.file), header = F)
sm <- fread(paste0("gunzip -cq ", sm.file))
gerp <- fread(paste0("gunzip -cq ", gerp.file))
cm <- fread(paste0("gunzip -cq ", cm.file))
ann <- fread(paste0("gunzip -cq ", ann.file))

### FORMAT

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

# continuous reference
colnames(ref) <- c("REF")
ref$POS <- 1:nrow(ref)

# calculate kmers
b1 <- ref$REF[1:(length(ref$REF)-6)]
b2 <- ref$REF[2:(length(ref$REF)-5)]
b3 <- ref$REF[3:(length(ref$REF)-4)]
b4 <- ref$REF[4:(length(ref$REF)-3)]
b5 <- ref$REF[5:(length(ref$REF)-2)]
b6 <- ref$REF[6:(length(ref$REF)-1)]
b7 <- ref$REF[7:(length(ref$REF))]
k7mer <- paste0(b1, b2, b3, b4, b5, b6, b7)
# k7mer <- stri_c(b1, b2, b3, b4, b5, b6, b7, sep = "_")
k7mer[grep("N", k7mer)] <- NA
ref$ANCESTRAL_7MER <- c(rep(NA, 3), k7mer, rep(NA, 3)) # ensure output equal length to input
rm(b1, b2, b3, b4, b5, b6, b7, k7mer)

### kmer counts

# for each annotation:
# subset positions and remove masked
# count 7mers

# subset masked int 
mask_int <- unique(unlist(seq2(from = masked$start, to = masked$end)))
k7_mu_all <- kmer_changes(7) # Identify all possible 11-mer substitutions

out_list <- list()
for (i in 1:length(unique(annotations$category))){
  
  # subset unmasked annotation kmers
  ann_sub <- annotations[category == unique(annotations$category)[i]] # subset annotation
  ann_int <- unique(unlist(seq2(from = ann_sub$start, to = ann_sub$end))) # subset annotation ind that are unmasked
  ann_int <- ann_int[!ann_int %in% mask_int]  # remove mask int from annotation int
  ref_sub <- ref[ann_int,]
  ref_sub <- ref_sub[complete.cases(ref_sub),]
  rm(ann_sub, ann_int)
  
  # count ancestral kmers
  k7_all <- data.table(k7_from = unique(k7_mu_all$k7_from))
  k7mer_from <- as.data.table(table(ref_sub$ANCESTRAL_7MER))
  colnames(k7mer_from) <- c(colnames(k7_all), "k7_from_N")
  k7mer_from <- k7mer_from[k7_all, on = "k7_from"]
  k7mer_from$k7_from_N[is.na(k7mer_from$k7_from_N)] <- 0
  k7mer_from$category <- unique(annotations$category)[i]

  out_list[[i]] <- k7mer_from
  print(unique(annotations$category)[i])  
}

output_all <- do.call("rbind", out_list)
output <- complement(output_all)
output <- output[c(which(stri_sub(output$k7_from, 4, -4) == "A"),
                   which(stri_sub(output$k7_from, 4, -4) == "C")), ]
output$chromosome <- chr
output <- output[,c("chromosome", "category", "k7_from", "k7_from_N")]

### OUTPUT
fwrite(x = output, file = out.file, append = T, compress = "gzip")

#######
