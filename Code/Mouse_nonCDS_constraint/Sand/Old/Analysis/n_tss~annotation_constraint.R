rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

# Correlation between regulatory feature constraint (UTRs, promoters, enhancers, CTCF binding sites) 
# and number of proximal genes.


### FUNCTION
# ann_cs_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_CS.csv"
# target_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_target_canTSSall.csv"
# canonical_file = "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_mouse_canPC_pos.csv"
# interest = c("Promoter", "Enhancer - proximal", "Enhancer - distal", "CTCF binding")
alakazam <- function(ann_cs_file,
                     target_file,
                     canonical_file,
                     interest){
  
  ### IMPORT
  ann_cs <- fread(ann_cs_file)
  target <- fread(target_file)
  canonical <- fread(canonical_file)
  
  ### FORMAT
  
  # subset columns
  target <- target[,c("chromosome", "category_start", "category_end", "category", "category_ID", "ensembl_transcript_id")]
  ann_cs <- ann_cs[,c("chromosome", "start", "end", "annotation", "ID", "length", "n_SNV", "n_sm", "n_CG", 
                      "exp_SNV", "residual", "OE_ratio", "OE_ratio_rank")]
  
  # colnames
  colnames(target) <- c("chromosome", "start", "end", "annotation", "ID", "ensembl_transcript_id")
  
  # subset annotations
  ann_cs <- ann_cs[annotation %in% interest]
  target <- target[annotation %in% interest]
  
  # subset canonical transcripts
  canonical <- unique(canonical$ensembl_transcript_id)
  target <- target[ensembl_transcript_id %in% canonical]
  
  # identify annotations with 0 transcripts
  no_TSS <- unique(ann_cs$ID[which(!ann_cs$ID %in% target$ID)])
  
  # calculate constraint percentiles for each annotation
  tmp_list <- list()
  for (cat in 1:length(unique(ann_cs$annotation))){
    
    tmp_sub <- ann_cs[annotation == unique(ann_cs$annotation)[cat]]
    
    percentile <- ecdf(tmp_sub$residual[!duplicated(tmp_sub$ID)])
    tmp_sub$residual_percentile <- percentile(tmp_sub$residual)
    
    percentile <- ecdf(tmp_sub$OE_ratio[!duplicated(tmp_sub$ID)])
    tmp_sub$OE_ratio_percentile <- percentile(tmp_sub$OE_ratio)
    
    percentile <- ecdf(tmp_sub$OE_ratio_rank[!duplicated(tmp_sub$ID)])
    tmp_sub$OE_ratio_rank_percentile <- percentile(tmp_sub$OE_ratio_rank)
    
    tmp_list[[cat]] <- tmp_sub
  }
  ann_cs <- do.call("rbind", tmp_list)
  rm(tmp_list, tmp_sub)
  
  # bin percentiles
  ann_cs$residual_percentile <- ceiling(ann_cs$residual_percentile * 100)
  ann_cs$OE_ratio_percentile <- ceiling(ann_cs$OE_ratio_percentile * 100)
  ann_cs$OE_ratio_rank_percentile <- ceiling(ann_cs$OE_ratio_rank_percentile * 100)
  
  # n TSS by annotation (including annotations with no TSS)
  dt <- target[ann_cs, on = c("chromosome", "start", "end", "annotation", "ID")]   # merge on ann_cs
  dt <- as.data.table(dt[complete.cases(dt),])   # remove NA
  n_tss <- dt[,c("annotation", "ID", "ensembl_transcript_id")]
  n_tss <- unique(n_tss)
  n_tss_tmp <- as.data.table(table(n_tss$ID))
  colnames(n_tss_tmp) <- c("ID", "n_tss")
  no_tss_tmp <- ann_cs[ID %in% no_TSS]
  no_tss_tmp <- no_tss_tmp[,c("annotation", "ID")]
  no_tss_tmp <- unique(no_tss_tmp)
  no_tss_tmp$n_tss <- 0
  n_tss <- n_tss[,c("annotation", "ID")]
  n_tss <- n_tss[n_tss_tmp, on = c("ID")]
  n_tss <- unique(n_tss)
  n_tss <- rbind(n_tss, no_tss_tmp)
  rm(n_tss_tmp, no_tss_tmp)
  
  # merge n_tss and annotation constraint
  dt_plot <- merge(n_tss, ann_cs[,c("ID", "residual_percentile", "OE_ratio_percentile", "OE_ratio_rank_percentile")])

  return(dt_plot)
}

### RUN 

m_plot <- alakazam(ann_cs_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_CS.csv",
                   target_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_target_canTSSall.csv",
                   canonical_file = "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_mouse_canPC_pos.csv",
                   interest = c("Promoter", "Enhancer - proximal", "Enhancer - distal", "CTCF binding"))

h_plot <- alakazam(ann_cs_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_CS.csv",
                   target_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_target_canTSSall.csv",
                   canonical_file = "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_pos.csv",
                   interest = c("Promoter", "Enhancer - proximal", "Enhancer - distal", "CTCF binding"))

fwrite(m_plot, "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_n_tss~annotation_constraint.csv")
fwrite(h_plot, "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Human_n_tss~annotation_constraint.csv")

#####






