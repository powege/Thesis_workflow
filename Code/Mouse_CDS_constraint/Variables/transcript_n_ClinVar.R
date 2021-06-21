# Script that counts synonymous, missense, nonsense, and spliceAD SNVs for each canonical transcript in QCed vcf

# SNV annotation categories -- 
# (https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)
# Non-functional:
# synonymous annotatiosn = "synonymous_variant", "stop_retained_variant","start_retained_variant"
# Functional: 
# missense annotations = "missense_variant"
# nonsnese annotations = "stop_gained", "start_lost", "stop_lost"
# splice annotations = "splice_donor_variant", "splice_acceptor_variant"

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

### SET PATHS

consequence_file_path <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP_canPC_IMPACT.vcf.gz"
out_file_path <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/ClinVar_germline_pathogenic_benign_canPC_n_SNV.csv"

### FUNCTIONS

# Function that counts the total synonymous SNVs
synonymous_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "synonymous_variant" |
                               names(consequence.vec) %like% "stop_retained_variant" |
                               names(consequence.vec) %like% "start_retained_variant"])
  return(out)
}

# Function that counts the total missense SNVs
missense_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "missense_variant" |
                               names(consequence.vec) %like% "protein_altering_variant"])
  return(out)
}

# Function that counts the total nonsense SNVs
nonsense_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "stop_gained" |
                               names(consequence.vec) %like% "start_lost" |
                               names(consequence.vec) %like% "stop_lost"])
  return(out)
}

# Function that counts the total spliceAD SNVs
spliceAD_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "splice_donor_variant" |
                               names(consequence.vec) %like% "splice_acceptor_variant"])
  return(out)
}

# Function that counts the total SNVs
total_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[ names(consequence.vec) %like% "synonymous_variant" |
                              names(consequence.vec) %like% "stop_retained_variant" |
                              names(consequence.vec) %like% "start_retained_variant" |
                              names(consequence.vec) %like% "missense_variant" |
                              names(consequence.vec) %like% "protein_altering_variant" |
                              names(consequence.vec) %like% "stop_gained" |
                              names(consequence.vec) %like% "start_lost" |
                              names(consequence.vec) %like% "stop_lost"])
                              
  return(out)
}

# Function that counts the non-functional SNVs, functional SNVs, and intron SNVs per Ensembl transcript
SNV_count <- function(data, colname){
  
  ### Subset HGNC_name, Ensembl_gene_id, Ensembl_transcript_id by unique Ensembl_transcript_id 
  gene_detail <- data[, c("Feature","SYMBOL")]
  gene_detail <- gene_detail[!duplicated(gene_detail$Feature),]
  colnames(gene_detail) <- c("ensembl_transcript_id", "external_gene_name")
  
  # ### Count synonymous SNVs
  # synonymous_all <- ddply(data, "Feature", synonymous_count)
  # colnames(synonymous_all) <- c("ensembl_transcript_id", "n_synonymous")
  # 
  # ### Count missense SNVs 
  # missense_all <- ddply(data, "Feature", missense_count)
  # colnames(missense_all) <- c("ensembl_transcript_id", "n_missense")
  # 
  # ### Count nonsense SNVs 
  # nonsense_all <- ddply(data, "Feature", nonsense_count)
  # colnames(nonsense_all) <- c("ensembl_transcript_id", "n_nonsense")
  # 
  # ### Count spliceAD SNVs 
  # spliceAD_all <- ddply(data, "Feature", spliceAD_count)
  # colnames(spliceAD_all) <- c("ensembl_transcript_id", "n_spliceAD")
  
  ### Count total SNVs 
  total_all <- ddply(data, "Feature", total_count)
  colnames(total_all) <- c("ensembl_transcript_id", colname)
  
  ### combine outputs
  output <- gene_detail
  # output <- merge(output, synonymous_all, all = T)
  # output <- merge(output, missense_all, all = T)
  # output <- merge(output, nonsense_all, all = T)
  # output <- merge(output, spliceAD_all, all = T)
  output <- merge(output, total_all, all = T)
  
  
  return(output)
}

### IMPORT
m.con <- fread(paste0("gunzip -cq ", consequence_file_path))

### FORMAT

# add column names
colnames(m.con) <- c("chromosome_name", "POS", "ID", "REF", "ALT", "status", "germline",
                      "Gene", "Feature", "Feature_type", "Consequence", "IMPACT", 
                      "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS")

# subset germline 
m.con <- m.con[germline %in% c("germline", "germline/somatic")]

### Get nonfunctional, functional and intron SNV counts (any MAF) per Ensembl transcript ID
pathogenic_counts <- SNV_count(m.con[status %like% "Pathogenic"], colname = "n_pathogenic")
benign_counts <- SNV_count(m.con[status %like% "Benign"], colname = "n_benign")
SNV_counts <- merge(pathogenic_counts, benign_counts, all = T, by = c("ensembl_transcript_id", "external_gene_name"))
SNV_counts$n_pathogenic[is.na(SNV_counts$n_pathogenic)] <- 0
SNV_counts$n_benign[is.na(SNV_counts$n_benign)] <- 0

# Get total genes
total_ens_tran <- length(unique(SNV_counts$ensembl_transcript_id))
total_HGNC <- length(unique(SNV_counts$external_gene_name))

# Get total SNV counts
# total_synonymous <- sum(SNV_counts$n_synonymous)
# total_missense <- sum(SNV_counts$n_missense)
# total_nonsense <- sum(SNV_counts$n_nonsense)
# total_spliceAD <- sum(SNV_counts$n_spliceAD)
total_pathogenic <- sum(SNV_counts$n_pathogenic)
total_benign <- sum(SNV_counts$n_benign)

SNV_summary <- data.frame(n_transcripts = total_ens_tran,
                            n_benign = total_benign,
                            n_pathogenic = total_pathogenic)

### OUTPUT
fwrite(SNV_counts, out_file_path)
# fwrite(SNV_summary, summary_file_path)



#####

# rm(list = ls())
# 
# library(data.table)
# library(tidyr)
# 
# dt <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps.vcf")
# colnames(dt) <- c("chromosome", "position", "id", "reference", "alternate", "significance", "germline", "type", "external_gene_name")
# dt <- separate_rows(dt, external_gene_name, convert = T, sep = ";")
# dt_pathogenic <- as.data.table(table(dt$external_gene_name[dt$significance %like% "Pathogenic" |
#                             dt$significance %like% "pathogenic"  ]))
# dt_benign <- as.data.table(table(dt$external_gene_name[dt$significance %like% "Benign" |
#                                                              dt$significance %like% "benign"  ]))
# colnames(dt_pathogenic) <- c("external_gene_name", "n_pathogenic")
# colnames(dt_benign) <- c("external_gene_name", "n_benign")
# SNV_counts <- merge(dt_pathogenic, dt_benign, all = T, by = c("external_gene_name"))
# fwrite(SNV_counts, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/ClinVar_germline_pathogenic_benign_n_snps.csv")

#####

# consequence_file_path <- "~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_PC_mouse_ortholog_o2o.csv"
# out_file_path <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/ClinVar_germline_pathogenic_benign_canPC_n_SNV.csv"
# m.con <- fread(consequence_file_path)
# pathogenic_counts <- SNV_count(m.con[Pathogenicity %like% "Pathogenic" |
#                                        Pathogenicity %like% "pathogenic" ], colname = "n_pathogenic")
# benign_counts <- SNV_count(m.con[Pathogenicity %like% "Benign" |
#                                    Pathogenicity %like% "benign"], colname = "n_benign")
# SNV_counts <- merge(pathogenic_counts, benign_counts, all = T, by = c("ensembl_transcript_id", "external_gene_name"))
# SNV_counts$n_pathogenic[is.na(SNV_counts$n_pathogenic)] <- 0
# SNV_counts$n_benign[is.na(SNV_counts$n_benign)] <- 0
# fwrite(SNV_counts, out_file_path)

