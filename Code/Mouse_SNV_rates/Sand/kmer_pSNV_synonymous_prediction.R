rm(list = ls())
graphics.off()

library(data.table)
library(robustbase)

################
### IMPORT DATA 
###############

n_snv <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/WM_Harr_etal_2016_allSPECIES_snps_PASS_VEP_canPC_n_SNV.csv")
k7_psnv <- fread("~/Dropbox/WM_Harr_etal_2016_allSPECIES_transcript_k7_pSNV.csv")
k5_psnv <- fread("~/Dropbox/WM_Harr_etal_2016_allSPECIES_transcript_k5_pSNV.csv")

##########
### FORMAT
##########

colnames(k7_psnv) <- c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "cds_length",           
                       "k7_p_synonymous", "k7_p_missense", "k7_p_nonsense")
colnames(k5_psnv) <- c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "cds_length",           
                       "k5_p_synonymous", "k5_p_missense", "k5_p_nonsense")

df <- k5_psnv[k7_psnv, on = c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "cds_length")]
df <- n_snv[df, on = c("external_gene_name", "ensembl_transcript_id")]
df <- unique(df[complete.cases(df),])

df$k7_p_nonsynonymous <- df$k7_p_missense + df$k7_p_nonsense
df$k5_p_nonsynonymous <- df$k5_p_missense + df$k5_p_nonsense

########################
### CALCULATE CONSTRAINT
#######################

## Fit model for synonymous variants
k7_mod <- lm(n_synonymous ~ k7_p_synonymous, data = df) # simple linear regression (ordinary least squares)
k7_coef <- summary(k7_mod)$coef
k7_r2 <- summary(k7_mod)$adj.r.squared

k5_mod <- lm(n_synonymous ~ k5_p_synonymous, data = df) # simple linear regression (ordinary least squares)
k5_coef <- summary(k5_mod)$coef
k5_r2 <- summary(k5_mod)$adj.r.squared

