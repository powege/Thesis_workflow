# SCRIPT that calculates constraint scores (O/E, Z-score, dN/dS)

rm(list = ls())
graphics.off()

library(data.table)
library(robustbase)

################
### IMPORT DATA 
###############

n_snv <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/WM_Harr_etal_2016_allSPECIES_snps_PASS_VEP_canPC_n_SNV.csv")
p_snv <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/WM_Harr_etal_2016_allSPECIES_transcript_k7_pSNV.csv")
coverage_m <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/WM_Harr_etal_2016_allSPECIES_transcript_n_coverage_masked.csv")
soft_m <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/m.mus_grc38_ensembl_transcript_n_soft_masked.csv")
meth <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/m.mus_grc38_ensembl_transcript_methylation.csv")
musc_nor <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_r.nor_orthologues.csv")
musc_spret <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_m.spret_orthologues.csv")


##########
### FORMAT
##########

meth$percent_CGunmethylated <- (meth$n_CGunmethylated / meth$n_CG) * 100

df <- n_snv[p_snv, on = c("external_gene_name", "ensembl_transcript_id")]
df <- coverage_m[df, on = c("ensembl_transcript_id")]
df <- soft_m[df, on = c("ensembl_transcript_id")]
df <- meth[df, on = c("ensembl_transcript_id")]
df <- unique(df[complete.cases(df),])

df$percent_soft_mask <- (df$n_soft_mask / df$cds_length) * 100
df$percent_dp_mask <- (df$n_dp_mask / df$cds_length) * 100
df$n_nonsynonymous <- df$n_missense + df$n_nonsense
df$p_nonsynonymous <- df$p_missense + df$p_nonsense

########################
### CALCULATE CONSTRAINT
#######################

## Set variables
y <- as.numeric(df$n_synonymous)
X_syn <- as.matrix(df[,c("p_synonymous",
                         # "percent_dp_mask",
                         "n_dp_mask",
                         # "percent_CGunmethylated"
                         "n_CGunmethylated"
                         )])
X_nonsyn <- as.matrix(df[,c("p_nonsynonymous",
                            # "percent_dp_mask",
                            "n_dp_mask",
                            # "percent_CGunmethylated"
                            "n_CGunmethylated"
                            )])

## Fit model for synonymous variants

# simple 
psyn_mod <- lm(n_synonymous ~ p_synonymous, data = df) # simple linear regression (ordinary least squares)
psyn_coef <- summary(psyn_mod)$coef
psyn_r2 <- summary(psyn_mod)$adj.r.squared

dp_mod <- lm(n_synonymous ~ n_dp_mask, data = df) # simple linear regression (ordinary least squares)
dp_coef <- summary(dp_mod)$coef
dp_r2 <- summary(dp_mod)$adj.r.squared

cg_mod <- lm(n_synonymous ~ percent_CGunmethylated, data = df) # simple linear regression (ordinary least squares)
cg_coef <- summary(cg_mod)$coef
cg_r2 <- summary(cg_mod)$adj.r.squared

mod2 <-  lm(n_synonymous ~ p_synonymous + n_dp_mask, data = df)
mod2_coef <- summary(mod2)$coef
mod2_r2 <- summary(mod2)$adj.r.squared

# multiple
syn_mod <- lm(y ~ X_syn) # multiple linear regression (ordinary least squares)
# plot(syn_mod)
coef <- summary(syn_mod)$coef
r2 <- summary(syn_mod)$adj.r.squared
# syn_mod <- lmrob(y ~ X_syn) # multiple linear regression (MM estimation to account for heteroskedastisity)
# plot(syn_mod_rob)
# coef <- summary(syn_mod)$coef
# r2 <- summary(syn_mod)$adj.r.squared

## Predict
syn_pred <- as.data.table(predict(syn_mod, 
                                  data.frame(X_syn), 
                                  interval="predict", 
                                  level=0.95))
colnames(syn_pred) <- c("exp_synonymous", "exp_synonymous_lwr", "exp_synonymous_upr")
df <- cbind(df, syn_pred)
df$exp_synonymous_lwr[which(df$exp_synonymous_lwr < 0)] <- 0

nonsyn_pred <- as.data.table(predict(syn_mod, 
                                      data.frame(X_nonsyn), 
                                      interval="predict", 
                                      level=0.95))
colnames(nonsyn_pred) <- c("exp_nonsynonymous", "exp_nonsynonymous_lwr", "exp_nonsynonymous_upr")
df <- cbind(df, nonsyn_pred)
df$exp_nonsynonymous_lwr[which(df$exp_nonsynonymous_lwr < 0)] <- 0

#################
### CALCULATE OE
################

df$exp_synonymous[df$exp_synonymous <= 0] <- 0
df$exp_synonymous_lwr[df$exp_synonymous_lwr <= 0] <- 0
df$exp_synonymous_upr[df$exp_synonymous_upr <= 0] <- 0
df$OE_synonymous <- df$n_synonymous/df$exp_synonymous
df$OE_synonymous_lwr <- df$n_synonymous/df$exp_synonymous_lwr
df$OE_synonymous_upr <- df$n_synonymous/df$exp_synonymous_upr
df$OE_synonymous[which(df$n_synonymous == 0 & df$exp_synonymous == 0)] <- 1
df$OE_synonymous_lwr[which(df$n_synonymous == 0 & df$exp_synonymous_lwr == 0)] <- 1
df$OE_synonymous_upr[which(df$n_synonymous == 0 & df$exp_synonymous_upr == 0)] <- 1

df$exp_nonsynonymous[df$exp_nonsynonymous <= 0] <- 0
df$exp_nonsynonymous_lwr[df$exp_nonsynonymous_lwr <= 0] <- 0
df$exp_nonsynonymous_upr[df$exp_nonsynonymous_upr <= 0] <- 0
df$OE_nonsynonymous <- df$n_nonsynonymous/df$exp_nonsynonymous
df$OE_nonsynonymous_lwr <- df$n_nonsynonymous/df$exp_nonsynonymous_lwr
df$OE_nonsynonymous_upr <- df$n_nonsynonymous/df$exp_nonsynonymous_upr
df$OE_nonsynonymous[which(df$n_nonsynonymous == 0 & df$exp_nonsynonymous <= 0)] <- 1
df$OE_nonsynonymous_lwr[which(df$n_nonsynonymous == 0 & df$exp_nonsynonymous_lwr == 0)] <- 1
df$OE_nonsynonymous_upr[which(df$n_nonsynonymous == 0 & df$exp_nonsynonymous_upr == 0)] <- 1

#####################
### CALCULATE z-score
####################

dif_syn <- df$exp_synonymous - df$n_synonymous
dif_syn_lwr <- df$exp_synonymous_lwr - df$n_synonymous
dif_syn_upr <- df$exp_synonymous_upr - df$n_synonymous
df$Z_synonymous <- (dif_syn - mean(dif_syn, na.rm = T))/sd(dif_syn, na.rm = T)
df$Z_synonymous_lwr <- (dif_syn_lwr - mean(dif_syn, na.rm = T))/sd(dif_syn, na.rm = T)
df$Z_synonymous_upr <- (dif_syn_upr - mean(dif_syn, na.rm = T))/sd(dif_syn, na.rm = T)

dif_nonsyn <- df$exp_nonsynonymous - df$n_nonsynonymous
dif_nonsyn_lwr <- df$exp_nonsynonymous_lwr - df$n_nonsynonymous
dif_nonsyn_upr <- df$exp_nonsynonymous_upr - df$n_nonsynonymous
df$Z_nonsynonymous <- (dif_nonsyn - mean(dif_nonsyn, na.rm = T))/sd(dif_nonsyn, na.rm = T)
df$Z_nonsynonymous_lwr <- (dif_nonsyn_lwr - mean(dif_nonsyn, na.rm = T))/sd(dif_nonsyn, na.rm = T)
df$Z_nonsynonymous_upr <- (dif_nonsyn_upr - mean(dif_nonsyn, na.rm = T))/sd(dif_nonsyn, na.rm = T)

####################
### CALCULATE dN/dS
##################

musc_nor <- musc_nor[orthology_type == "ortholog_one2one"]
musc_nor$rat_dNdS <- musc_nor$dn/musc_nor$ds
musc_nor$rat_dNdS[musc_nor$dn == 0 & musc_nor$ds == 0] <- 1
musc_nor <- musc_nor[,c("M_ensembl_transcript_id", "rat_dNdS")]
colnames(musc_nor) <- c("ensembl_transcript_id", "rat_dNdS")

musc_spret <- musc_spret[orthology_type == "ortholog_one2one"]
musc_spret$spretus_dNdS <- musc_spret$dn/musc_spret$ds
musc_spret$rat_dNdS[musc_spret$dn == 0 & musc_spret$ds == 0] <- 1
musc_spret <- musc_spret[,c("musculus_ensembl_transcript_id", "spretus_dNdS")]
colnames(musc_spret) <- c("ensembl_transcript_id", "spretus_dNdS")

df <- musc_nor[df, on = "ensembl_transcript_id"]
df <- musc_spret[df, on = "ensembl_transcript_id"]


###########
### OUTPUT
##########

output <- df[,c(
                "external_gene_name",
                "ensembl_transcript_id",
                "cds_length",
                "n_CG",
                "n_CGunmethylated",       
                "n_synonymous",
                "n_nonsynonymous",          
                "p_synonymous",
                "p_nonsynonymous",
                "percent_soft_mask",
                "percent_dp_mask",      
                "exp_synonymous",
                # "exp_synonymous_lwr",
                # "exp_synonymous_upr",
                "exp_nonsynonymous",
                # "exp_nonsynonymous_lwr",
                # "exp_nonsynonymous_upr",
                # "OE_synonymous_lwr",
                "OE_synonymous", 
                # "OE_synonymous_upr",
                # "OE_nonsynonymous_lwr",
                "OE_nonsynonymous",
                # "OE_nonsynonymous_upr",
                # "Z_synonymous_lwr",
                "Z_synonymous", 
                "Z_synonymous_upr",
                # "Z_nonsynonymous_lwr",
                "Z_nonsynonymous",
                # "Z_nonsynonymous_upr",
                "rat_dNdS",
                "spretus_dNdS"
                )]
fwrite(output, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Constraint_scores/WM_Harr_etal_2016_allSPECIES_constraint_scores.csv")




#####
