rm(list = ls())
graphics.off()

library(data.table)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

in.file <- args[1]
out.file <- args[2]

### IMPORT
vcf <- fread(in.file)


### SPLIT INFO COLUMN
vcf[, c("Gene","Feature","Feature_type","Consequence","IMPACT","SYMBOL","SYMBOL_SOURCE","BIOTYPE","CANONICAL","CCDS") := tstrsplit(INFO, "|", fixed=TRUE)][, c("INFO"):=NULL]

### REMOVE CSQ=
vcf$Gene <- gsub("CSQ=", "", vcf$Gene)

### SUBSET PROTEIN CODING CANONICAL 
vcf <- vcf[BIOTYPE == "protein_coding" & CANONICAL == "YES"]

### SUBSET BY IMPACT
vcf <- vcf[IMPACT == "LOW" | IMPACT == "MODERATE" | IMPACT == "HIGH"]

### EXPORT
fwrite(vcf, 
       out.file,
       compress = "gzip")



