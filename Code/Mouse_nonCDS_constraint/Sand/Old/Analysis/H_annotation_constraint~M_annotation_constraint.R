rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

# Correlation in constraint between regulatory features (UTRs, promoters, enhancers, CTCF binding sites)
# of orthologous genes. 

### SET ARGS 
interest <- c("CTCF binding",
              "Enhancer - distal",
              "Enhancer - proximal", 
              "Promoter", 
              "Exon - UTR")

### IMPORT
m_ann_cs <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_CS.csv")
h_ann_cs <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_CS.csv")

m_target <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_target_canTSSone.csv")
h_target <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_target_canTSSone.csv")

m_utr <- fread("~/Dropbox/PhD/Data/Target_gene/Mouse_UTR_tranID.csv")
h_utr <- fread("~/Dropbox/PhD/Data/Target_gene/Human_UTR_tranID.csv")

m_ann_all <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv")
h_ann_all <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv")

orths <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_HM_orthologues.csv")

### FORMAT

# subset columns
m_utr <- m_utr[,c("V1",  "V6", "V7", "V5", "V2")]
h_utr <- h_utr[,c("V1",  "V6", "V7", "V5", "V2")]

m_ann_all <- m_ann_all[,c("V1", "V2", "V3", "V4", "V6")]
h_ann_all <- h_ann_all[,c("V1", "V2", "V3", "V4", "V6")]

m_target <- m_target[,c("chromosome", "category_start", "category_end", "category", "category_ID", "ensembl_transcript_id")]
h_target <- h_target[,c("chromosome", "category_start", "category_end", "category", "category_ID", "ensembl_transcript_id")]

m_ann_cs <- m_ann_cs[,c("chromosome", "start", "end", "annotation", "ID", "length", "n_SNV", "n_sm", "n_CG", 
                    "exp_SNV", "residual", "OE_ratio")]
h_ann_cs <- h_ann_cs[,c("chromosome", "start", "end", "annotation", "ID", "length", "n_SNV", "n_sm", "n_CG", 
                    "exp_SNV", "residual", "OE_ratio")]

orths <- orths[,c("H_external_gene_name", "M_external_gene_name", "H_ensembl_transcript_id",
                  "M_ensembl_transcript_id", "orthology_type")]

# colnames
colnames(m_utr) <- paste0("M_", c("chromosome", "start", "end", "annotation", "ensembl_transcript_id"))
colnames(h_utr) <- paste0("H_", c("chromosome", "start", "end", "annotation", "ensembl_transcript_id"))

colnames(m_ann_all) <- paste0("M_", c("chromosome", "start", "end", "annotation", "ID"))
colnames(h_ann_all) <- paste0("H_", c("chromosome", "start", "end", "annotation", "ID"))

colnames(m_target) <- paste0("M_", c("chromosome", "start", "end", "annotation", "ID", "ensembl_transcript_id"))
colnames(h_target) <- paste0("H_", c("chromosome", "start", "end", "annotation", "ID", "ensembl_transcript_id"))

colnames(m_ann_cs) <- paste0("M_", colnames(m_ann_cs))
colnames(h_ann_cs) <- paste0("H_", colnames(h_ann_cs))

# merge
m_dt <- m_ann_all[m_utr, on = paste0("M_", c("chromosome", "start", "end", "annotation"))]
h_dt <- h_ann_all[h_utr, on = paste0("H_", c("chromosome", "start", "end", "annotation"))]

m_dt <- rbind.fill(m_target, m_dt)
h_dt <- rbind.fill(h_target, h_dt)

m_dt <- merge(m_dt, m_ann_cs, all = T)
h_dt <- merge(h_dt, h_ann_cs, all = T)

m_dt <- merge(m_dt, orths, by = paste0("M_", "ensembl_transcript_id"))
h_dt <- merge(h_dt, orths, by = paste0("H_", "ensembl_transcript_id"))

dt <- merge(m_dt, h_dt, all = F)
rm(m_ann_cs, m_target, m_utr, m_ann_all, h_ann_cs, h_target, h_utr, h_ann_all, orths, h_dt, m_dt)

# QC
dt <- as.data.table(dt[complete.cases(dt),])
dt <- unique(dt)

# calculate constraint percentiles for each mouse annotation
tmp_list <- list()
for (cat in 1:length(unique(dt$M_annotation))){
  tmp_sub <- dt[M_annotation == unique(dt$M_annotation)[cat]]
  
  percentile <- ecdf(tmp_sub$M_residual[!duplicated(tmp_sub$M_ID)])
  tmp_sub$M_residual_percentile <- percentile(tmp_sub$M_residual)
  
  percentile <- ecdf(tmp_sub$M_OE_ratio[!duplicated(tmp_sub$M_ID)])
  tmp_sub$M_OE_ratio_percentile <- percentile(tmp_sub$M_OE_ratio)
  
  tmp_list[[cat]] <- tmp_sub
}
dt <- do.call("rbind", tmp_list)
rm(tmp_list, tmp_sub)

# calculate constraint percentiles for each human annotation
tmp_list <- list()
for (cat in 1:length(unique(dt$H_annotation))){
  tmp_sub <- dt[H_annotation == unique(dt$H_annotation)[cat]]
  
  percentile <- ecdf(tmp_sub$H_residual[!duplicated(tmp_sub$H_ID)])
  tmp_sub$H_residual_percentile <- percentile(tmp_sub$H_residual)
  
  percentile <- ecdf(tmp_sub$H_OE_ratio[!duplicated(tmp_sub$H_ID)])
  tmp_sub$H_OE_ratio_percentile <- percentile(tmp_sub$H_OE_ratio)
  
  tmp_list[[cat]] <- tmp_sub
}
dt <- do.call("rbind", tmp_list)
rm(tmp_list, tmp_sub)

# bin annotation constraint percentiles
dt$M_residual_percentile <- ceiling(dt$M_residual_percentile * 100)
dt$M_OE_ratio_percentile <- ceiling(dt$M_OE_ratio_percentile * 100)

dt$H_residual_percentile <- ceiling(dt$H_residual_percentile * 100)
dt$H_OE_ratio_percentile <- ceiling(dt$H_OE_ratio_percentile * 100)

# subset annotations
dt_plot <- dt[M_annotation %in% interest]
dt_plot <- dt_plot[M_annotation == H_annotation]
dt_plot <- unique(dt_plot)

### EXPORT 
fwrite(dt_plot, "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/H_annotation_constraint~M_annotation_constraint.csv")


#####






