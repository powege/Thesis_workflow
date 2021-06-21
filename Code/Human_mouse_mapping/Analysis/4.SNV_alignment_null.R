rm(list=ls())
graphics.off()

library(data.table)

# INPUT:
# human annotation bed 
# mouse annotation bed
# human mouse alignment bed
# human mendelian vcf
# human gwas vcf
# chromosome

# OUTPUT:
# append bp counts for alignment and conservation

### SET VARS 

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments
if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} 

# set args variables
h.ann.file <- args[1]
m.ann.file <- args[2]
align.file <- args[3]
in.snv.file <- args[4]
chr <- as.character(args[5])
functions.file <- args[6]
out.null.file <- args[7]

# h.ann.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/h.sap_GRC38_v101_whole_genome_features_multicell.csv.gz"
# m.ann.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/m.mus_GRC38_v101_whole_genome_features_multicell.csv.gz"
# align.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/hsap_grch38_v_mmus_grcm38_v101_alignment_chr15.bed.gz"
# # in.snv.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_SNV_alignment_by_chr.csv"
# in.snv.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_by_chr.csv"
# chr <- as.character(15)
# functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"
# # out.null.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_SNV_alignment_NULL_by_chr.csv"
# out.null.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_NUL_by_chr.csv"

### FUNCTIONS
source(functions.file)

### IMPORT 
human_ann <- fread(paste0("gunzip -cq ", h.ann.file))
mouse_ann <- fread(paste0("gunzip -cq ", m.ann.file))
align <- fread(paste0("gunzip -cq ", align.file))
snv_obs <- fread(in.snv.file)

### FORMAT 

# subset by human chromosome
human_ann <- human_ann[chromosome == chr]
align <- align[chromosome_human == chr]
snv_obs <- snv_obs[chromosome == chr]

# QC
human_ann <- human_ann[percentage_Nmask < 1]
mouse_ann <- mouse_ann[percentage_Nmask < 1]

### RUN FOR ALL ANNNOTATION

# snv_bed = snv_obs
alakazam <- function(human_ann, mouse_ann, align, snv_bed, out_file){
# for each human annotation in chromosome:
ann <- unique(human_ann$category)
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  # sample snv_bed$n_total[snv_bed$annotation == ann[i]] pos from human_ann
  snv_sample <- bed.sample.pos(bed = human_ann[category == ann[i]][,c("chromosome", "start", "end")], nPOS = snv_bed$n_total[snv_bed$annotation == ann[i]])
  n_total[i] <- nrow(unique(snv_sample))
  
  if (n_total[i] != 0){
  # orthologous sequences for human annotation to get total alignment
  snv_ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                    "chromosome_mouse", "start_mouse", "end_mouse")],
                               bed = snv_sample)
  n_aligned[i] <- nrow(unique(snv_ann_align[,c("chromosome_A", "start_A", "end_A")]))
  
  # total of mouse annotation in orthologous sequences for conservation
  snv_ann_conserved <- unique(bed.intersect(bed1 = snv_ann_align[,c("chromosome_B", "start_B", "end_B")],
                                   bed2 = mouse_ann[category == ann[i]][,c("chromosome", "start", "end")]))
  n_conserved[i] <- nrow(snv_ann_conserved)
  } else {
    n_aligned[i] <- 0
    n_conserved[i] <- 0
  }
  
  print(ann[i])
}

output <- data.table( chromosome = chr,
                      annotation = ann,
                      n_total = n_total,
                      n_aligned = n_aligned,
                      n_conserved = n_conserved)
fwrite(output, out_file, append = T)
# rm(snv_ann, snv_ann_align, snv_ann_conserved, ann, i, n_aligned, n_conserved, n_total)
}

replicate(200, alakazam(human_ann, mouse_ann, align, snv_obs, out.null.file), simplify = FALSE)


### RUN FOR "Functional - all" "Non-CDS" and "Total"

# snv_bed <- snv_obs
kadabra <- function(human_ann, mouse_ann, align, snv_bed, out_file){
human_ann_2_list <- list(functional = human_ann[category %in% c("Exon - CDS",
                                                                "Exon - UTR",
                                                                "Exon - other",
                                                                "Promoter",           
                                                                "Enhancer - proximal",
                                                                "Enhancer - distal",
                                                                "CTCF binding",
                                                                "Miscellaneous",      
                                                                "TAD boundry",
                                                                "Intron - proximal")],
                         non_cds =  human_ann[!category %in% c("Exon - CDS", "Intron - distal", "Unannotated")],
                         total = human_ann)

mouse_ann_2_list <- list(functional = mouse_ann[category %in% c("Exon - CDS",
                                                                "Exon - UTR",
                                                                "Exon - other",
                                                                "Promoter",           
                                                                "Enhancer - proximal",
                                                                "Enhancer - distal",
                                                                "CTCF binding",
                                                                "Miscellaneous",      
                                                                "TAD boundry",
                                                                "Intron - proximal")],
                         non_cds =  mouse_ann[!category %in% c("Exon - CDS", "Intron - distal", "Unannotated")],
                         total = mouse_ann)

ann <- c("Functional - all", "Functional - nonCDS", "Total")
ann_list <- list(c("Exon - CDS",
                   "Exon - UTR",
                   "Exon - other",
                   "Promoter",           
                   "Enhancer - proximal",
                   "Enhancer - distal",
                   "CTCF binding",
                   "Miscellaneous",      
                   "TAD boundry",
                   "Intron - proximal"),
                 c("Exon - UTR",
                   "Exon - other",
                   "Promoter",           
                   "Enhancer - proximal",
                   "Enhancer - distal",
                   "CTCF binding",
                   "Miscellaneous",      
                   "TAD boundry",
                   "Intron - proximal"),
                 c("Exon - CDS",
                   "Exon - UTR",
                   "Exon - other",
                   "Promoter",           
                   "Enhancer - proximal",
                   "Enhancer - distal",
                   "CTCF binding",
                   "Miscellaneous",      
                   "TAD boundry",
                   "Intron - proximal",
                   "Intron - distal", 
                   "Unannotated"))
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  human_ann_2 <- human_ann_2_list[[i]]
  mouse_ann_2 <- mouse_ann_2_list[[i]]
  
  snv_sample <- bed.sample.pos(bed = human_ann_2[,c("chromosome", "start", "end")], nPOS = sum(snv_bed$n_total[snv_bed$annotation %in% ann_list[[i]]]))
  n_total[i] <- unique(nrow(snv_sample))
  
  if (n_total[i] != 0){
    # orthologous sequences for human annotation to get total alignment
    snv_ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                          "chromosome_mouse", "start_mouse", "end_mouse")],
                                     bed = snv_sample)
    n_aligned[i] <- nrow(unique(snv_ann_align[,c("chromosome_A", "start_A", "end_A")]))
    
    # total of mouse annotation in orthologous sequences for conservation
    # snv_ann_conserved <- unique(bed.intersect(bed1 = snv_ann_align[,c("chromosome_B", "start_B", "end_B")],
    #                                           bed2 = mouse_ann_2[,c("chromosome", "start", "end")]))
    # n_conserved[i] <- nrow(snv_ann_conserved)
  } else {
    n_aligned[i] <- 0
    # n_conserved[i] <- 0
  }
  print(ann[i])
}
output <- data.table( chromosome = chr,
                        annotation = ann,
                        n_total = n_total,
                        n_aligned = n_aligned,
                        n_conserved = n_conserved)

fwrite(output, out_file, append = T)
}

replicate(200, kadabra(human_ann, mouse_ann, align, snv_obs, out.null.file), simplify = FALSE)




