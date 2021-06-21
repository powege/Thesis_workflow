### SCRIPT that calculates the total of each sequence in a bed file that matches a given vector
### INPUT: BED file of sequences; vector of POS
### OUTPUT: BED file of sequences with: n_sm (soft mask); 
                                      # n_Nm (N mask); 
                                      # n_cm (coverage mask);

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
sm_file <- args[3]
Nm_file <- args[4]
cm_file <- args[5]
chr <- args[6]

# bed_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- ""
# sm_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_smPOS_Ensembl_GRC38_v94_chr19.csv"
# Nm_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_NPOS_Ensembl_GRC38_v94_chr19.csv"
# cm_file <- "~/Dropbox/PhD/Data/MGP/bam/MGP_fraction90_depth_less10X_chr19.txt"
# chr <- "19"

# bed_file <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# sm_file <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_smPOS_Ensembl_GRC38_v94_chr21.csv"
# Nm_file <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_NPOS_Ensembl_GRC38_v94_chr21.csv"
# cm_file <- "/well/lindgren/George/Data/gnomAD/coverage/gnomad_v2.1_fraction90_depth_less10X_chr21.tsv"
# out_file <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_sNc_masking_chr21.csv"
# chr <- "19"


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
sm_pos <- fread(sm_file, header = F)
Nm_pos <- fread(Nm_file, header = F)
cm_pos <- fread(cm_file, header = F)

### FORMAT

# subset chomosome
bed <- bed[V1 == chr]
# remove incomplete cases
bed <- subset(bed, !is.na(bed$V2) & !is.na(bed$V3))
# remove duplicated rows
bed <- bed[!duplicated(bed),]
# subset start and end
bed_sub <- bed[,c(2,3)]

# convert to vectors
cm_pos <- as.integer(unique(cm_pos$V1))
Nm_pos <- as.integer(unique(Nm_pos$V1))
sm_pos <- as.integer(unique(sm_pos$V1))

# apply function
bed$n_sm <- bed_match(bed = bed_sub, vec = sm_pos)
bed$n_Nm <- bed_match(bed = bed_sub, vec = Nm_pos)
bed$n_cm <- bed_match(bed = bed_sub, vec = cm_pos)

### EXPORT
fwrite(bed, out_file, col.names = F)


#####

### STACK OVERFLOW



