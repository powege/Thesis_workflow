rm(list = ls())
graphics.off()

library(data.table)

### FUNNCTIONS
# dt_in <- dt[,c("mmus_external_gene_name",
#                "OE_nonsynonymous_percentile",
#                "n_pathogenic_kb",
#                "n_benign_kb"
#                )]
dugtrio <- function(dt_in){
  
  colnames(dt_in) <- c("mmus_external_gene_name",
                       "score",
                       "n_pathogenic_kb",
                       "n_benign_kb"
                       )
  dt_in <- dt_in[complete.cases(dt_in)]
  percentile <- ecdf(dt_in$score[!duplicated(dt_in$mmus_external_gene_name)])
  dt_in$score <- ceiling(percentile(dt_in$score)*100)
  
  # plot(dt_in$n_pathogenic_kb ~ dt_in$score)
  # mod <- glm(formula = n_pathogenic_kb ~ score, data = dt_in, family = quasibinomial("logit"))
  # summary(mod)
  
  # dt_ratio <- aggregate(list(average_ratio=dt_in$ratio), by=list(score_percentile=dt_in$score), FUN=mean)
  # plot(dt_ratio$average_ratio~dt_ratio$score_percentile)
  # cor.test(dt_ratio$average_ratio, dt_ratio$score_percentile)
  
  dt_pathogenic <- aggregate(list(average_pathogenic=dt_in$n_pathogenic_kb), by=list(score_percentile=dt_in$score), FUN=mean)
  plot(dt_pathogenic$average_pathogenic~dt_pathogenic$score_percentile)
  cor.test(dt_pathogenic$average_pathogenic, dt_pathogenic$score_percentile)
  
  dt_benign <- aggregate(list(average_benign=dt_in$n_benign_kb), by=list(score_percentile=dt_in$score), FUN=mean)
  plot(dt_benign$average_benign~dt_benign$score_percentile)
  cor.test(dt_benign$average_benign, dt_benign$score_percentile)
  
  return(list(dt_in, dt_mn, dt_md))
}

### SET VARS
infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv"
# outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/"

### IMPORT
dt <- fread(infile)

### FORMAT
dt <- dt[,c("mmus_external_gene_name",
            "Z_nonsynonymous",
            "OE_nonsynonymous",
            "rat_dNdS",
            "spretus_dNdS",
            "pLI",
            "oe_lof_upper",
            "mis_z",
            "Z_nonsynonymous_percentile",
            "OE_nonsynonymous_percentile",
            "rat_dNdS_percentile",
            "spretus_dNdS_percentile",
            "pLI_percentile",
            "oe_lof_upper_percentile",
            "mis_z_percentile",
            "cds_length",
            "n_pathogenic",               
            "n_benign",
            "hsap_orthology_type")]
dt <- dt[!is.na(dt$cds_length)]
dt <- dt[hsap_orthology_type == "ortholog_one2one"]
dt$n_pathogenic[is.na(dt$n_pathogenic)] <- 0
dt$n_benign[is.na(dt$n_benign)] <- 0
# dt <- dt[!is.na(dt$n_pathogenic) | !is.na(dt$n_benign)]
dt$n_pathogenic_kb <- (dt$n_pathogenic / dt$cds_length) * 1000 
dt$n_benign_kb <- (dt$n_benign / dt$cds_length) * 1000 
# dt$n_pathogenic_benign_ratio <- dt$n_pathogenic / dt$n_benign
# dt$ratio <- (((dt$n_pathogenic +1) / (dt$n_benign +1)) / dt$cds_length) *1000



#####

# infile <- "~/Dropbox/PhD/Data/PC_constraint/Supplementary_data_table.csv"
# df <- fread(infile)
# x <- df[M_external_gene_name == "Brca1"]
# df <- df[,c("M_external_gene_name",
#             "Z_nonsynonymous",
#             "H_cds_length",
#             "n_pathogenic",
#             "n_benign")]
# df <- df[!is.na(df$H_cds_length)]
# df$n_pathogenic_kb <- (df$n_pathogenic / df$H_cds_length) * 1000 
# df$n_benign_kb <- (df$n_benign / df$H_cds_length) * 1000 
# 
# dt_in <- df[,c("M_external_gene_name",
#                "Z_nonsynonymous",
#                "n_pathogenic_kb",
#                "n_benign_kb"
#                )]
# 
# 
# infile <- "~/Dropbox/PhD/Data/PC_constraint/Supplementary_data_table.csv"
# df <- fread(infile)
# df <- df[,c("M_external_gene_name",
#             "n_pathogenic",
#             "n_benign")]
# colnames(df) <- c("mmus_external_gene_name", 
#                   "n_pathogenic2",
#                   "n_benign2")
# infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv"
# dt <- fread(infile)
# dt <- dt[,c("mmus_external_gene_name",
#             "n_pathogenic",
#             "n_benign")][df, on = "mmus_external_gene_name"]
# plot(dt$n_pathogenic, dt$n_pathogenic2)
# plot(dt$n_benign, dt$n_benign2)


# rm(list=ls())
# 
# cv <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_benign_snps_QCed_VEP_v94.vcf", fill = T)
# cv2 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP_canPC_IMPACT.vcf.gz")
# 
# 
# head(cv)
# head(cv2)
# 
# source("~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R")
# 
# cv$end <- cv$V2
# cv2$end <- cv2$V2
# x <- bed.intersect(cv2[,c("V1", "V2", "end")], cv[,c("V1", "V2", "end")])
# 
# length(unique(cv2$V3))
# table(cv2$V3 %in% cv$V3)
