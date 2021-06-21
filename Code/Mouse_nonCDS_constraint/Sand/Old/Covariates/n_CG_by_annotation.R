### SCRIPT that calculates the number of CG dinucleotides given a BED file

rm(list = ls())
graphics.off()

library(data.table)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) { stop("arguments must be supplied", call.=FALSE) } 

### SET ARGS

bed_file <- args[1]
ref_file <- args[2]
snv_file <- args[3]
out_file <- args[4]
chr <- as.integer(args[5])

# bed_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv"
# ref_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr21.csv"
# SNV_file <- "~/Dropbox/PhD/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_controls_chr21.vcf"
# out_file <- ""
# chr <- 21

### FUNCTIONS

# vectorised version of substr
substr2 <- Vectorize(substr, vectorize.args = c("start", "stop"))

### IMPORT

ref_dt <- fread(ref_file)
cat <- fread(cat_file)


### FORMAT

colnames(cat) <- c("chromosome", "start", "end", "category", "strand", "ID") # colnames
cat <- cat[chromosome == chr] # subset chromosome
# cat <- cat[category %in% interest] # subset annotations of interest
cat <- unique(cat) # remove duplicates

colnames(ref_dt) <- c("POS", "REF")

if (ref_dt$POS[1] == 1 & ref_dt$POS[nrow(ref_dt)] == nrow(ref_dt)){ # if reference genome continuous
ref <- paste0(ref_dt$REF, collapse = "") # collapse referecne to string
}

cat$sequence <- unlist(substr2(x = ref, start = cat$start, stop = cat$end)) # extract sequence
cat$n_CG <- (str_count(cat$sequence, "CG") + str_count(cat$sequence, "GC"))
cat$sequence <- NULL


### EXPORT
fwrite(cat, out_file, col.names = F)


#####
