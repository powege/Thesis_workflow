# SCRIPT that creates data table of variables by region for each chromosome (regions defined by distance between SNV sites)

rm(list = ls())
graphics.off()

library(data.table)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {  stop("Arguments must be supplied", call.=FALSE) } 

# set args
species <- args[1]
chromo <- as.integer(args[2])
ref.file <- args[3]
vcf.file <- args[4]
smask.file <- args[5]
Nmask.file <- args[6]
coverage.file <- args[7]
out.file <- args[8]

# species <- "human"
# chromo <- 21
# ref.file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr21.csv"
# vcf.file <- "~/Dropbox/PhD/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_controls_chr21.vcf"
# smask.file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_smPOS_Ensembl_GRC38_v94_chr21.csv"
# Nmask.file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_NPOS_Ensembl_GRC38_v94_chr21.csv"
# coverage.file <- "~/Dropbox/PhD/Data/gnomAD/coverage/gnomad_v2.1_fraction90_depth_less10X_chr21.tsv"
# out.file <- "~/Dropbox/PhD/Data/NC_constraint/Variables/Human_constrainted_regions_gnomAD_variables_chr21.csv"

### FUNCTIONS

# FUNCTION vectorised version of substr
substr2 <- Vectorize(substr, vectorize.args = c("start", "stop"))

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

# reference genome
ref <- fread(ref.file, fill = T)
# SNVs
vcf <- fread(vcf.file, fill = T)
# soft mask
smask <- fread(smask.file, fill = T)
# N mask
Nmask <- fread(Nmask.file, fill = T)
# coverage
coverage <- fread(coverage.file, fill = T)


### FORMAT

# colnames
colnames(ref) <- c("POS", "REF")

# check reference genome continuous
if (ref$POS[1] == 1 & ref$POS[nrow(ref)] == nrow(ref)){ 
  print("reference continuous") 
} else {
  print("reference  not continuous") 
}

# subset MAF
if (species == "human"){ vcf <- vcf[vcf$V19 >= 0.001,] }  

# identify POS with SNV
vcf <- as.integer(unique(vcf$V2))
vcf <- vcf[!is.na(vcf)]

# identify POS with N
Nmask <- as.integer(unique(Nmask$V1))
Nmask <- Nmask[!is.na(Nmask)]

# combine vcf and Nmask for break POS
breaks <- sort(unique(c(vcf, Nmask)))

# convert vector to bed file 
bed <- data.table(chromosome = rep(chromo, (length(breaks) + 1)),
                  start = c(1, (breaks + 1)),
                  end = c((breaks - 1), nrow(ref)))
# bed <- data.table(chromosome = rep(chromo, ((length(breaks) * 2) + 1)),
#                   start = c(1, (breaks + 1), breaks),
#                   end = c((breaks - 1), nrow(ref), breaks))
# bed <- bed[order(start),]

# remove rows where start > end 
bed <- bed[start < end]

# collapse referecne to string
ref <- paste0(ref$REF, collapse = "") 

# fraction CG
bed$sequence <- unlist(substr2(x = ref, start = bed$start, stop = bed$end)) 
bed$n_CG <- (str_count(bed$sequence, "CG") + str_count(bed$sequence, "GC"))
bed$sequence <- NULL
bed$f_CG <- bed$n_CG / ((bed$end + 1) - bed$start)

# subset bed 
bed_sub <- bed[,c(2,3)]

# # fraction with Nmask
# Nmask <- as.integer(unique(Nmask$V1))
# Nmask <- Nmask[!is.na(Nmask)]
# bed$n_Nm <- bed_match(bed = bed_sub, vec = Nmask)
# bed$f_Nm <- bed$n_Nm / ((bed$end + 1) - bed$start)

# fraction with soft mask
smask <- as.integer(unique(smask$V1))
smask <- smask[!is.na(smask)]
bed$n_sm <- bed_match(bed = bed_sub, vec = smask)
bed$f_sm <- bed$n_sm / ((bed$end + 1) - bed$start)

# fraction with low coverage
coverage <- as.integer(unique(coverage$V1))
coverage <- coverage[!is.na(coverage)]
bed$n_cm <- bed_match(bed = bed_sub, vec = coverage)
bed$f_cm <- bed$n_cm / ((bed$end + 1) - bed$start)

# coverage weighted length
bed$length <- (bed$end + 1) - bed$start
bed$length_cm_weighted <- bed$length * (1 - bed$f_cm)

output <- bed[,c("chromosome", "start", "end", "f_CG", "f_sm", "f_cm", "length_cm_weighted")]

### EXPORT  
fwrite(output, out.file, col.names = T)


#####

### STACK OVRFLOW

# library(data.table)
# set.seed(1)
# vec <- sample(1:10000, 1000)
# vec <- sort(vec)
# bed <- data.table(start = c(1, vec),
#                   end = c((vec - 1), 10000))
# bed <- data.table(start = c(1, (vec + 1), vec),
#                   end = c((vec - 1), 10000, vec))
# bed <- bed[order(start),]


