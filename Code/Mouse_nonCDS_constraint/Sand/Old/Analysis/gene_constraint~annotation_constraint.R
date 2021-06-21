rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(ggplot2)

# Correlation in constraint between regulatory features (UTRs, promoters, enhancers, CTCF binding sites) 
# and closest gene (oe ratio).

### RUN FOR MOUSE

### SET ARGS 
prefix <- "M_"
interest <- c("CTCF binding",
              "Enhancer - distal",
              "Enhancer - proximal", 
              "Promoter", 
              "Exon - UTR")

### IMPORT
ann_cs <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_CS.csv")
target <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_target_canTSSone.csv")
utr <- fread("~/Dropbox/PhD/Data/Target_gene/Mouse_UTR_tranID.csv")
ann_all <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv")
gene_cs <- fread("~/Dropbox/PhD/Data/Constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt")
orths <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_HM_orthologues.csv")

### FORMAT

# subset columns
utr <- utr[,c("V1",  "V6", "V7", "V5", "V2")]
ann_all <- ann_all[,c("V1", "V2", "V3", "V4", "V6")]
target <- target[,c("chromosome", "category_start", "category_end", "category", "category_ID", "ensembl_transcript_id")]
ann_cs <- ann_cs[,c("chromosome", "start", "end", "annotation", "ID", "length", "n_SNV", "n_sm", "n_CG", 
                    "exp_SNV", "residual", "OE_ratio")]
gene_cs <- gene_cs[,c("gene", "transcript", "pLI", "oe_lof_upper", "lof_z")]
orths <- orths[,c("H_external_gene_name", "M_external_gene_name", "H_ensembl_transcript_id",
                  "M_ensembl_transcript_id", "orthology_type")]

# colnames
colnames(utr) <- paste0(prefix, c("chromosome", "start", "end", "annotation", "ensembl_transcript_id"))
colnames(ann_all) <- paste0(prefix, c("chromosome", "start", "end", "annotation", "ID"))
colnames(target) <- paste0(prefix, c("chromosome", "start", "end", "annotation", "ID", "ensembl_transcript_id"))
colnames(ann_cs) <- paste0(prefix, colnames(ann_cs))
colnames(gene_cs) <- paste0("H_", c("external_gene_name", "ensembl_transcript_id", "pLI", "oe_lof_upper", "lof_z"))

# merge
dt <- ann_all[utr, on = paste0(prefix, c("chromosome", "start", "end", "annotation"))]
dt <- rbind.fill(target, dt)
dt <- merge(dt, ann_cs, all = T)
dt <- merge(dt, orths, by = paste0(prefix, "ensembl_transcript_id"))
dt <- merge(dt, gene_cs, all = F)
rm(ann_cs, gene_cs, target, utr, orths, ann_all)

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

# calculate constraint percentiles for each huamn gene
percentile <- ecdf(dt$H_pLI[!duplicated(dt$H_ensembl_transcript_id)])
dt$H_pLI_percentile <- percentile(dt$H_pLI)
percentile <- ecdf(dt$H_oe_lof_upper[!duplicated(dt$H_ensembl_transcript_id)])
dt$H_oe_lof_upper_percentile <- percentile(dt$H_oe_lof_upper)
percentile <- ecdf(dt$H_lof_z[!duplicated(dt$H_ensembl_transcript_id)])
dt$H_lof_z_percentile <- percentile(dt$H_lof_z)

# bin annotation constraint percentiles
dt$M_residual_percentile <- ceiling(dt$M_residual_percentile * 100)
dt$M_OE_ratio_percentile <- ceiling(dt$M_OE_ratio_percentile * 100)

# subset annotations
dt_plot <- dt[M_annotation %in% interest]

### EXPORT 
fwrite(dt_plot, "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_gene_constraint~annotation_constraint.csv")


### RUN FOR HUMAN

### SET ARGS 
prefix <- "H_"
interest <- c("CTCF binding",
              "Enhancer - distal",
              "Enhancer - proximal", 
              "Promoter", 
              "Exon - UTR")

### IMPORT
ann_cs <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_CS.csv")
target <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_target_canTSSone.csv")
utr <- fread("~/Dropbox/PhD/Data/Target_gene/Human_UTR_tranID.csv")
ann_all <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv")
gene_cs <- fread("~/Dropbox/PhD/Data/Constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt")
orths <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_HM_orthologues.csv")

### FORMAT

# subset columns
utr <- utr[,c("V1",  "V6", "V7", "V5", "V2")]
ann_all <- ann_all[,c("V1", "V2", "V3", "V4", "V6")]
target <- target[,c("chromosome", "category_start", "category_end", "category", "category_ID", "ensembl_transcript_id")]
ann_cs <- ann_cs[,c("chromosome", "start", "end", "annotation", "ID", "length", "n_SNV", "n_sm", "n_CG", 
                    "exp_SNV", "residual", "OE_ratio")]
gene_cs <- gene_cs[,c("gene", "transcript", "pLI", "oe_lof_upper", "lof_z")]
orths <- orths[,c("H_external_gene_name", "M_external_gene_name", "H_ensembl_transcript_id",
                  "M_ensembl_transcript_id", "orthology_type")]

# colnames
colnames(utr) <- paste0(prefix, c("chromosome", "start", "end", "annotation", "ensembl_transcript_id"))
colnames(ann_all) <- paste0(prefix, c("chromosome", "start", "end", "annotation", "ID"))
colnames(target) <- paste0(prefix, c("chromosome", "start", "end", "annotation", "ID", "ensembl_transcript_id"))
colnames(ann_cs) <- paste0(prefix, colnames(ann_cs))
colnames(gene_cs) <- paste0("H_", c("external_gene_name", "ensembl_transcript_id", "pLI", "oe_lof_upper", "lof_z"))

# merge
dt <- ann_all[utr, on = paste0(prefix, c("chromosome", "start", "end", "annotation"))]
dt <- rbind.fill(target, dt)
dt <- merge(dt, ann_cs, all = T)
dt <- merge(dt, orths, by = paste0(prefix, "ensembl_transcript_id"))
dt <- merge(dt, gene_cs, all = F)
rm(ann_cs, gene_cs, target, utr, orths, ann_all)

# QC
dt <- as.data.table(dt[complete.cases(dt),])
dt <- unique(dt)

# calculate constraint percentiles for each mouse annotation
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

# calculate constraint percentiles for each huamn gene
percentile <- ecdf(dt$H_pLI[!duplicated(dt$H_ensembl_transcript_id)])
dt$H_pLI_percentile <- percentile(dt$H_pLI)
percentile <- ecdf(dt$H_oe_lof_upper[!duplicated(dt$H_ensembl_transcript_id)])
dt$H_oe_lof_upper_percentile <- percentile(dt$H_oe_lof_upper)
percentile <- ecdf(dt$H_lof_z[!duplicated(dt$H_ensembl_transcript_id)])
dt$H_lof_z_percentile <- percentile(dt$H_lof_z)

# bin annotation constraint percentiles
dt$H_residual_percentile <- ceiling(dt$H_residual_percentile * 100)
dt$H_OE_ratio_percentile <- ceiling(dt$H_OE_ratio_percentile * 100)

# subset annotations
dt_plot <- dt[H_annotation %in% interest]

### EXPORT 
fwrite(dt_plot, "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Human_gene_constraint~annotation_constraint.csv")


#####







