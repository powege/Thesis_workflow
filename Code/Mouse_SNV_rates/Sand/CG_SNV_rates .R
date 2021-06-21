rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS
source("~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R")

### SET VARS
sm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
cm.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr"
gerp.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/gerp_constrained_elements.mus_musculus.bed.gz"
ann.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
meth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz"
snv.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_CG_snps_PASS.vcf.gz"

### IMPORT
sm <- fread(paste0("gunzip -cq ", sm.file))
gerp <- fread(paste0("gunzip -cq ", gerp.file))
cm <- list()
for (i in 1:19){
  cm[[i]] <- fread(paste0("gunzip -cq ", cm.file, i, ".bed.gz"))
}
cm <- do.call("rbind", cm)
ann <- fread(paste0("gunzip -cq ", ann.file))
cgs <- fread(paste0("gunzip -cq ", meth.file))
snvs <- fread(paste0("gunzip -cq ", snv.file))

# creadt mask bed
colnames(sm) <- c("chromosome", "start", "end")
colnames(gerp) <- c("chromosome", "start", "end")
colnames(cm) <- c("chromosome", "start", "end")
ann <- ann[category %in% c("Exon - CDS", "Exon - UTR", "Exon - other", "Promoter", "Enhancer - proximal", "Enhancer - distal", "CTCF binding", "Miscellaneous", "Intron - proximal")][,c("chromosome", "start", "end")]
mask <- collapse.overlap(rbind(sm, gerp, cm, ann))

cgs_back <- cgs
ind <- bed.intersect(cgs[,c("chromosome", "start", "end")], mask) # ind is the masked CGs
ind$chromosome <- as.integer(ind$chromosome)
cgs <- cgs[!ind, on = c("chromosome", "start", "end")] # antijoinn

cg_m <- cgs[coverage > 5 & percent_methylated > 60]
cg_um <- cgs[coverage > 5 & percent_methylated < 20]
# cg_m <- cgs[-which(cgs$coverage > 5 & cgs$percent_methylated < 20)]
# cg_um <- cgs[which(cgs$coverage > 5 & cgs$percent_methylated < 20)]

snv_CT <- rbind( snvs[V4 == "C" & V5 == "T"], snvs[V4 == "G" & V5 == "A"])
snv_CA <- rbind( snvs[V4 == "C" & V5 == "A"], snvs[V4 == "G" & V5 == "T"])
snv_CG <- rbind( snvs[V4 == "C" & V5 == "G"], snvs[V4 == "G" & V5 == "C"])

# all CG C > T
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CT[V1 == i]
  cg_sub <- cgs[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_all <- do.call("rbind", out_list)
CT_all <- nrow(snv_all) / (nrow(cgs)*2)

# meth CG C > T
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CT[V1 == i]
  cg_sub <- cg_m[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_m <- do.call("rbind", out_list)
CT_m <- nrow(snv_m) / (nrow(cg_m)*2)

# unmeth CG C > T
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CT[V1 == i]
  cg_sub <- cg_um[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_um <- do.call("rbind", out_list)
CT_um <- nrow(snv_um) / (nrow(cg_um)*2)
 
# all CG C > A
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CA[V1 == i]
  cg_sub <- cgs[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_all <- do.call("rbind", out_list)
CA_all <- nrow(snv_all) / (nrow(cgs)*2)

# meth CG C > A
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CA[V1 == i]
  cg_sub <- cg_m[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_m <- do.call("rbind", out_list)
CA_m <- nrow(snv_m) / (nrow(cg_m)*2)

# unmeth CG C > A
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CA[V1 == i]
  cg_sub <- cg_um[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_um <- do.call("rbind", out_list)
CA_um <- nrow(snv_um) / (nrow(cg_um)*2)

# all CG C > G
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CG[V1 == i]
  cg_sub <- cgs[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_all <- do.call("rbind", out_list)
CG_all <- nrow(snv_all) / (nrow(cgs)*2)

# meth CG C > G
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CG[V1 == i]
  cg_sub <- cg_m[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_m <- do.call("rbind", out_list)
CG_m <- nrow(snv_m) / (nrow(cg_m)*2)

# unmeth CG C > G
out_list <- list()
for(i in 1:19){
  snv_sub <- snv_CG[V1 == i]
  cg_sub <- cg_um[chromosome == i]
  out_list[[i]] <- snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))]
  print(i)
}
snv_um <- do.call("rbind", out_list)
CG_um <- nrow(snv_um) / (nrow(cg_um)*2)


# differnece cound be due to ancestral inferance
# test improvment relative to the k1CG and k7 (no meth data)





