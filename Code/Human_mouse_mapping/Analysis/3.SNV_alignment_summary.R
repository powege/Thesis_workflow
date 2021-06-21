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
h.mendelian.file <- args[4]
h.gwas.file <- args[5]
chr <- as.character(args[6])
functions.file <- args[7]
out.mendelian.file <- args[8]
out.gwas.file <- args[9]

# h.ann.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/h.sap_GRC38_v101_whole_genome_features_multicell.csv.gz"
# m.ann.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/m.mus_GRC38_v101_whole_genome_features_multicell.csv.gz"
# align.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/hsap_grch38_v_mmus_grcm38_v101_alignment_chr15.bed.gz"
# h.mendelian.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps.vcf"
# h.gwas.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/GWAS/Formatted/GWAS_catalog_DDC_QCed_2020_11_30.tsv"
# chr <- as.character(15)
# functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"
# out.mendelian.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_ClinVar_Pathogenic_SNV_alignment_by_chr.csv"
# out.gwas.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_by_chr.csv"

### FUNCTIONS
source(functions.file)

### IMPORT 
human_ann <- fread(paste0("gunzip -cq ", h.ann.file))
mouse_ann <- fread(paste0("gunzip -cq ", m.ann.file))
align <- fread(paste0("gunzip -cq ", align.file))
mendel <- fread(h.mendelian.file)
gwas <- fread(h.gwas.file)

### FORMAT 

colnames(mendel) <- c("chromosome", "start", "dbSNP", "Reference", "Alternate",
                      "ClinicalSignificance", "OriginSimple")
mendel <- mendel[OriginSimple == "germline" & ClinicalSignificance %like% "Pathogenic"]
mendel <- mendel[chromosome == chr]
mendel <- mendel[,c("chromosome", "start")]
mendel$end <- mendel$start

colnames(gwas) <- c("chromosome", "start", "Parent_term")
gwas <- gwas[chromosome == chr]
gwas <- gwas[,c("chromosome", "start")]
gwas$end <- gwas$start

# subset by human chromosome
human_ann <- human_ann[chromosome == chr]
align <- align[chromosome_human == chr]

# QC
human_ann <- human_ann[percentage_Nmask < 1]
mouse_ann <- mouse_ann[percentage_Nmask < 1]

### RUN FOR ALL ANNNOTATION

alakazam <- function(human_ann, mouse_ann, align, snv_bed){
# for each human annotation in chromosome:
ann <- unique(human_ann$category)
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  snv_ann <- unique(bed.intersect(bed1 = snv_bed, bed2 = human_ann[category == ann[i]][,c("chromosome", "start", "end")]))
  n_total[i] <- nrow(snv_ann)
  
  if (n_total[i] != 0){
  # orthologous sequences for human annotation to get total alignment
  snv_ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                    "chromosome_mouse", "start_mouse", "end_mouse")],
                               bed = snv_ann)
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

return(data.table( chromosome = chr,
                          annotation = ann,
                          n_total = n_total,
                          n_aligned = n_aligned,
                          n_conserved = n_conserved))
rm(snv_ann, snv_ann_align, snv_ann_conserved, ann, i, n_aligned, n_conserved, n_total)
}

mendel_out_1 <- alakazam(human_ann, mouse_ann, align, mendel)
gwas_out_1 <- alakazam(human_ann, mouse_ann, align, gwas)


### RUN FOR "Functional - all" "Non-CDS" and "Total"

# snv_bed <- mendel
kadabra <- function(human_ann, mouse_ann, align, snv_bed){
human_ann_2_list <- list(functional = human_ann[category %in% c("Exon - CDS",
                                                                "Exon - 5'UTR",
                                                                "Exon - 3'UTR",
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
                                                                "Exon - 5'UTR",
                                                                "Exon - 3'UTR",
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
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  human_ann_2 <- human_ann_2_list[[i]]
  mouse_ann_2 <- mouse_ann_2_list[[i]]
  
  snv_ann <- unique(bed.intersect(bed1 = snv_bed, bed2 = human_ann_2[,c("chromosome", "start", "end")]))
  n_total[i] <- nrow(snv_ann)
  
  if (n_total[i] != 0){
    # orthologous sequences for human annotation to get total alignment
    snv_ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                          "chromosome_mouse", "start_mouse", "end_mouse")],
                                     bed = snv_ann)
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
return(data.table( chromosome = chr,
                        annotation = ann,
                        n_total = n_total,
                        n_aligned = n_aligned,
                        n_conserved = n_conserved))
}

mendel_out_2 <- kadabra(human_ann, mouse_ann, align, mendel)
gwas_out_2 <- kadabra(human_ann, mouse_ann, align, gwas)

mendel_out <- rbind(mendel_out_1, mendel_out_2)
gwas_out <- rbind(gwas_out_1, gwas_out_2)

### EXPORT

fwrite(mendel_out, out.mendelian.file, append = T)
fwrite(gwas_out, out.gwas.file, append = T)


