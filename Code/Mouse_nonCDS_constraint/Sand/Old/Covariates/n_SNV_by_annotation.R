### SCRIPT that calculates the total of each sequence in a bed file that matches a given vector
### INPUT: BED file of sequences; vector of POS
### OUTPUT: BED file of sequences with: n_SNV

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) { stop("arguments must be supplied", call.=FALSE) } 

### SET ARGS
bed_file <- args[1]
out_file <- args[2]
SNV_file <- args[3]
chr <- args[4]
species <- args[5]

# bed_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv"
# out_file <- ""
# SNV_file <- "~/Dropbox/PhD/Data/MGP/Variants/vcf_QCed_VEP/MGP_v5_allMUSMUS_snps_QCed_VEP_v94_chr19.vcf"
# chr <- "19"
# species <- "mouse"
# 
# bed_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv"
# out_file <- ""
# SNV_file <- "~/Dropbox/PhD/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_controls_chr21.vcf"
# chr <- "21"
# species <- "human"

### FUNCTIONS

### FUNCTION that returns the seqence match between bed file and vector
# (ie for each sequence in bed file1, how many integers overlap with vector?)
bed_match <- function(bed, vec){
  
  bed_dt <- setDT(bed) # set as bed as data.table
  colnames(bed_dt) <- c("start", "end")
  bed_dt[, ind := .I] # add uniqe index to data.table
  setkey(bed_dt) # sets keys // order data by all columns
  
  vec_dt <- as.data.table(vec, key = 'vec') # convert to data.table
  vec_dt[, vec2 := vec] # dublicate column
  
  # Fast overlap join:
  ans1 = foverlaps(vec_dt, bed_dt, by.x = c('vec', 'vec2'), by.y = c('start', 'end'),
                   type = "within", nomatch = 0L)
  counts <- ans1[, .N, keyby = ind] # count by ind
  # merge to inital data
  bed_dt[, n_match := counts[bed_dt, on = .(ind), x.N]]
  setorder(bed_dt, ind) # reorder by ind to get inital order
  bed_dt[, ind := NULL] # deletes ind colum
  bed_dt[is.na(n_match), n_match := 0L] # NAs is 0 count
  
  return(bed_dt$n_match) 
}

### IMPORT
bed <- fread(bed_file, header = F)
snv_pos <- fread(SNV_file, header = F, fill = T)

### FORMAT

# subset chomosome
bed <- bed[V1 == chr]
# remove incomplete cases
bed <- subset(bed, !is.na(bed$V2) & !is.na(bed$V3))
# remove duplicated rows
bed <- unique(bed)
# subset start and end
bed_sub <- bed[,c(2,3)]

# convert to vector
if (species == "mouse"){ snv_pos <- as.integer(unique(snv_pos$V2)) }
if (species == "human"){ snv_pos <- as.integer(unique(snv_pos$V2[snv_pos$V19 >= 0.001])) }
# remove NA
snv_pos <- snv_pos[!is.na(snv_pos)]

# apply function
bed$n_SNV <- bed_match(bed = bed_sub, vec = snv_pos)

### EXPORT
fwrite(bed, out_file, col.names = F)


#####

### STACK OVERFLOW



