rm(list=ls())
graphics.off()

library(data.table)

#### SET PATHS
seq.in.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_canPC_seq.csv"
# rem.out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/"
QCed.out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_canPC_seq_QCpass.csv"

### IMPORT 
seq <- fread(seq.in.file)

### QC 

removed <- data.frame()

# remove seq not divisible by 3
ncod <- seq$cds_length/3
rm.id <- which(ncod%%1!=0)
if (length(rm.id) != 0){
  removed <- rbind(removed, seq[rm.id,])
  seq <- seq[-rm.id,]
}

# remove seq with nchar < 9
rm.id <- which(seq$cds_length<9)
if (length(rm.id) != 0){
  removed <- rbind(removed, seq[rm.id,])
  seq <- seq[-rm.id,]
}

# remove seq with nchar > 15000
rm.id <- which(seq$cds_length>15000)
if (length(rm.id) != 0){
  removed <- rbind(removed, seq[rm.id,])
  seq <- seq[-rm.id,]
}

# # remove seq with nchar < 300
# rm.id <- which(seq$cds_length<300)
# if (length(rm.id) != 0){
#   removed <- rbind(removed, seq[rm.id,])
#   seq <- seq[-rm.id,]
# }

# remove seq with characters other than ATCG 
rm.id <- grep('[^ATGC]', seq$coding)
if (length(rm.id) != 0){
  removed <- rbind(removed, seq[rm.id,])
  seq <- seq[-rm.id,]
}

# remove seq that do not start with start codon
c1 <- substring(seq$coding, 1, 3)
rm.id <- which(c1 != "ATG")
if (length(rm.id) != 0){
  removed <- rbind(removed, seq[rm.id,])
  seq <- seq[-rm.id,]
}

# remove seq that do not end with stop codon
cn <- substring(seq$coding, seq$cds_length-2, seq$cds_length)
rm.id <- which(cn != "TAG" & cn != "TAA" & cn != "TGA")
if (length(rm.id) != 0){
  removed <- rbind(removed, seq[rm.id,])
  seq <- seq[-rm.id,]
}

### REORDER
seq <- seq[order(seq$chromosome_name)]
removed <- removed[order(removed$chromosome_name)]

seq <- seq[complete.cases(seq),]
removed <- removed[complete.cases(removed),]


### OUTPUT 
# fwrite(removed, rem.out.file)
fwrite(seq, QCed.out.file)

