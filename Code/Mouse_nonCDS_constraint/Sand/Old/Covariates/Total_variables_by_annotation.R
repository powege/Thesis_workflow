### SCRIPT that counts the total CG dinucleotides, SNVs, soft mask, by chromosome and annotation

rm(list = ls())
graphics.off()

library(data.table)
library(stringr)
library(dplyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) { stop("arguments must be supplied", call.=FALSE) } 

### SET ARGS

ref_file <- args[1]
bed_file <- args[2]
SNV_file <- args[3]
out_file <- args[4]
chr <- as.integer(args[5])
species <- args[6]

# ref_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr21.csv"
# bed_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv"
# SNV_file <- "~/Dropbox/PhD/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_controls_chr21.vcf"
# out_file <- ""
# chr <- 21
# species <- "human"

# ref_file <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr22.csv" 
# bed_file <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv" 
# SNV_file <- "/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_controls_chr22.vcf" 
# out_file <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_SNV_CG_total_chr22.csv" 
# chr <- "22"
# species <- "human"


### FUNCTIONS

# vectorised version of substr
substr2 <- Vectorize(substr, vectorize.args = c("start", "stop"))

# FUNCTION that collapses overlapping sequences
collapse_overlap <- function(sub){
  
  require(data.table)
  require(dplyr)
  
  # ensure start <= end
  dt <- sub[ start > end, `:=`( start = end, end = start)]
  
  dt2 <- rbind(dt,dt)
  dt2[1:(.N/2), ID := 1]
  dt2[(.N/2 +1):.N, ID := 2]
  setkey(dt2, ID, start, end)
  
  squished <- dt2[,.(START_DT = start, 
                     END_DT = end, 
                     indx = c(0, cumsum(as.numeric(lead(start)) > cummax(as.numeric(end)))[-.N])),
                  keyby=ID
                  ][,.(start=min(START_DT), 
                       end = max(END_DT)),
                    by=c("ID","indx")
                    ]
  squished <- squished[, c("start", "end")]
  squished <- unique(squished)
  
  return(squished[,c("start", "end")])
}

# FUNCTION that returns the seqence match between bed file and vector
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

bed <- fread(bed_file)
ref <- fread(ref_file)
snv <- fread(SNV_file, header = F, fill = T)

### FORMAT 

bed <- bed[,c("V1", "V2", "V3", "V4")]
colnames(bed) <- c("chromosome", "start", "end", "category")
colnames(ref) <- c("POS", "REF")

# colapse seqeunces
bed_chr <- bed[chromosome == chr]
ann_list <- list()
for (ann in 1:length(unique(bed$category))){
    sub_ann <- bed_chr[category == unique(bed$category)[ann]]
    out_ann <- collapse_overlap(sub_ann)
    out_ann$category <- unique(bed$category)[ann]
    ann_list[[ann]] <- out_ann
}
out_chr <- do.call("rbind", ann_list)
out_chr$chromosome <- chr

# reference vector to string
if (ref$POS[1] == 1 & ref$POS[nrow(ref)] == nrow(ref)){ # if reference genome continuous
  ref <- paste0(ref$REF, collapse = "") # collapse referecne to string
}

# n_CG
out_chr$sequence <- unlist(substr2(x = ref, start = out_chr$start, stop = out_chr$end)) # extract sequence
out_chr$n_CG <- (str_count(out_chr$sequence, "CG") + str_count(out_chr$sequence, "GC"))
out_chr$sequence <- NULL

# n_SNV
out_chr <- subset(out_chr, !is.na(out_chr$start) & !is.na(out_chr$end)) # remove incomplete cases
out_chr <- unique(out_chr) # remove duplicated rows
out_chr_sub <- out_chr[,c("start", "end")] # subset start and end
# convert to vector
if (species == "mouse"){ snv <- as.integer(unique(snv$V2)) }
if (species == "human"){ snv <- as.integer(unique(snv$V2[snv$V19 >= 0.001])) }
snv <- snv[!is.na(snv)] # remove NA
out_chr$n_SNV <- bed_match(bed = out_chr_sub, vec = snv) # apply function

# total var by chromosome
length_chr <- nchar(ref)
n_CG_chr <- (str_count(ref, "CG") + str_count(ref, "GC"))
n_SNV_chr <- length(snv)

# total var by annotation
length_ann <- rep(NA, length(unique(out_chr$category)))
n_CG_ann <- rep(NA, length(unique(out_chr$category)))
n_SNV_ann <- rep(NA, length(unique(out_chr$category)))
for (i in 1:length(unique(out_chr$category))){
  sub_tmp <- out_chr[category == unique(out_chr$category)[i]]
  length_ann[i] <- sum((sub_tmp$end + 1) - sub_tmp$start)
  n_CG_ann[i] <- sum(sub_tmp$n_CG)
  n_SNV_ann[i] <- sum(sub_tmp$n_SNV)
}

output_total <- data.table(chromosome = chr,
                           annotation = "total",
                           length = length_chr,
                           n_CG = n_CG_chr,
                           n_SNV = n_SNV_chr)
output_ann <- data.table(chromosome = rep(chr, length(unique(out_chr$category))),
                         annotation = unique(out_chr$category),
                         length = length_ann,
                         n_CG = n_CG_ann,
                         n_SNV = n_SNV_ann)
output <- rbind(output_ann, output_total)                         

### EXPORT
fwrite(output, out_file, col.names = F)
                         
                         