rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
source("~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R")

# Load data ========== 

sm <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz")

dp_list <- list()
for (i in 1:19){
  dp_list[[i]] <- fread(paste0("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr", i, ".bed.gz"))
}
dp <- do.call("rbind", dp_list)

ann <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz")

len <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_summary.bed")

# Format ==========

colnames(sm) <- c("chromosome", "start", "end")
colnames(dp) <- c("chromosome", "start", "end")
ann <- ann[category %in% c("Exon - CDS", "Exon - UTR", "Exon - other", "Promoter", "Enhancer - proximal", "Enhancer - distal", "CTCF binding", "Miscellaneous", "TAD boundry", "Intron - proximal")][, c("chromosome", "start", "end")]
len <- len[, 1:3]
colnames(len) <- c("chromosome", "start", "end")

mask <- rbind(sm, dp, ann)
mask <- collapse.overlap(mask)
mask$length <- (mask$end + 1) - mask$start
mask_len <- mask[, mask_length := sum(.SD), by=.(chromosome), .SDcols="length"]
# mask_len <- mask[, sum(.SD), by=.(chromosome), .SDcols="length"]
mask_len <- unique(mask_len[,c("chromosome", "mask_length")])
mask_len$chromosome <- as.integer(mask_len$chromosome)

len$length <- (len$end + 1) - len$start
chr_len <- len[,c("chromosome", "length")]

dt <- chr_len[mask_len, on = "chromosome"]
sum(dt$mask_length) / sum(dt$length)

dt$unmasked_length <- dt$length - dt$mask_length
dt <- dt[order(unmasked_length),]

sum(dt$unmasked_length[1:4]) / sum(dt$unmasked_length)
sum(dt$unmasked_length[16:19]) / sum(dt$unmasked_length)


#####

x <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1_pSNV.csv.gz")
sum(x$k1_to_N[x$k1_to == "*"])



