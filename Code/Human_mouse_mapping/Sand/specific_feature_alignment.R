# Given human gene name:
# Identify coordinates for all promoters; enhancers; and UTRs for human gene and mouse orthologue
# Calculate alignmnet and conservation of gene-specific features

rm(list=ls())
graphics.off()

library(data.table)
library(plyr)
library(matrixStats)

### FUNCTIONS

source("~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R")

# tss <- can[,c("transcription_start_site", "strand")]
# distal_extension = 300000
dominion <- function(tss, distal_extension){
  
  colnames(tss) <- c("ensembl_transcript_id", "transcription_start_site", "strand")
  
  # subset by strand 
  strand1 <- subset(tss, tss$strand == 1)
  strand2 <- subset(tss, tss$strand == -1)
  
  # get start and end pos 
  strand1$proximal_start <- strand1$transcription_start_site - 5000
  strand1$proximal_end <- strand1$transcription_start_site + 1000
  strand1$distal1_start <- strand1$transcription_start_site - (distal_extension + 5000)
  strand1$distal1_end <- strand1$transcription_start_site - 5001
  strand1$distal2_start <- strand1$transcription_start_site + 1001
  strand1$distal2_end <- strand1$transcription_start_site + (distal_extension + 1000)
  
  strand2$proximal_start <- strand2$transcription_start_site + 5000
  strand2$proximal_end <- strand2$transcription_start_site - 1000
  strand2$distal1_start <- strand2$transcription_start_site + (distal_extension + 5000)
  strand2$distal1_end <- strand2$transcription_start_site + 5001
  strand2$distal2_start <- strand2$transcription_start_site - 1001
  strand2$distal2_end <- strand2$transcription_start_site - (distal_extension + 1000)
  
  # rbind
  tss <- rbind(strand1, strand2)
  rm(strand1, strand2)
  
  # melt 
  id_vars <- c("ensembl_transcript_id","transcription_start_site","strand")
  dt_start <- tss[,c(id_vars,"proximal_start", "distal1_start","distal2_start"), with = F]
  dt_start <- melt(dt_start, id.vars = id_vars)
  colnames(dt_start) <- c(id_vars, "domain", "domain_start")
  dt_start$domain <- as.character(dt_start$domain)
  dt_start$domain[dt_start$domain == "proximal_start"] <- "proximal"
  dt_start$domain[dt_start$domain == "distal1_start"] <- "distal1"
  dt_start$domain[dt_start$domain == "distal2_start"] <- "distal2"
  dt_start <- unique(dt_start)
  
  dt_end <- tss[,c(id_vars,"proximal_end", "distal1_end","distal2_end"), with = F]
  dt_end <- melt(dt_end, id.vars = id_vars)
  colnames(dt_end) <- c(id_vars, "domain", "domain_end")
  dt_end$domain <- as.character(dt_end$domain)
  dt_end$domain[dt_end$domain == "proximal_end"] <- "proximal"
  dt_end$domain[dt_end$domain == "distal1_end"] <- "distal1"
  dt_end$domain[dt_end$domain == "distal2_end"] <- "distal2"
  dt_end <- unique(dt_end)
  
  tss <- dt_start[dt_end, on = c(id_vars, "domain")]
  rm(dt_start, dt_end)
  
  return(tss)
}


# reg <- ann[,c("chromosome",
#               "start",
#               "end",
#               "category_id")]
# tss <- can[, c("chromosome_name",
#                "ensembl_transcript_id",
#                "transcription_start_site",
#                "domain",
#                "domain_start",
#                "domain_end" )]
# chr <- 1
all_canonical <- function(reg, tss, chr){
  
  colnames(reg) <- c("chromosome",
                     "category_start",
                     "category_end",
                     "category_id",
                     "regulatory_feature_stable_id")
  colnames(tss) <- c("chromosome",
                     "ensembl_transcript_id",
                     "transcription_start_site",
                     "domain",
                     "domain_start",
                     "domain_end" )
  
  # subset chromosome
  reg_sub <- subset(reg, reg$chromosome == chr)
  tss_sub <- subset(tss, tss$chromosome == chr)
  
  # ensure start <= end
  reg_sub[ category_start > category_end, `:=`( category_start = category_end, category_end = category_start)]
  tss_sub[ domain_start > domain_end, `:=`( domain_start = domain_end, domain_end = domain_start)]
  
  # merge by overlap
  setkey(reg_sub, chromosome, category_start, category_end)
  setkey(tss_sub, chromosome, domain_start, domain_end)
  ans <- unique( foverlaps( tss_sub, reg_sub ))
  
  # calculate total overlap
  ans[, overlap_start := rowMaxs( as.matrix(.SD), na.rm = TRUE ), .SDcols = c("category_start", "domain_start")]
  ans[, overlap_end   := rowMins( as.matrix(.SD), na.rm = TRUE ), .SDcols = c("category_end", "domain_end")]
  ans[, overlap_size  := overlap_end - overlap_start + 1 ]
  
  # complete cases
  ans <- ans[complete.cases(ans),]
  
  # calculate min distance to TSS from cat start and end
  ans$TSS_dist <- apply( data.table( start_dist = abs( (ans$transcription_start_site + 1) - ans$category_start ),
                                     end_dist = abs( (ans$transcription_start_site + 1) - ans$category_end) ), 
                         1,
                         min
  )
  
  # overlap fractions
  ans$category_length <- (ans$category_end + 1) - ans$category_start
  ans$domain_length <- (ans$domain_end + 1) - ans$domain_start
  ans$f_category_overlap <- ans$overlap_size/ans$category_length
  ans$f_domain_overlap <- ans$overlap_size/ans$domain_length
  
  # subset variables
  ans <- ans[,c("chromosome", 
                "category_start", 
                "category_end", 
                "category_id", 
                "regulatory_feature_stable_id",
                "ensembl_transcript_id", 
                "transcription_start_site",
                "TSS_dist", 
                "domain", 
                "f_category_overlap", 
                "f_domain_overlap")]
  ans <- ans[order(category_id, category_start, category_end),] # order 
  # ans <- ans[ensembl_transcript_id %in% can_ID] # subset canonical
  
  return(ans)
}

# FUNCTION that returns row with shortest TSS_dist
# sub <- output_all[category_ID == "D1"]
closest_one <- function(sub){ sub[which(sub$TSS_dist == min(sub$TSS_dist)),] }

# FUNCTION that prioritises proximal overlap
prioritise_proximal <- function(sub){ 
  if ("proximal" %in% sub$domain) {
    out <- sub[which(sub$domain == "proximal"),]
  } else ( out <- sub )
  return(out)
}



onix <- function(can, ann, gene_name){
  
# subset transcript of interest
can <- can[external_gene_name == gene_name]

can <- unique(can[,c("chromosome_name",
                     "external_gene_name",
                     "ensembl_gene_id",
                     "ensembl_transcript_id",
                     "transcription_start_site",
                     "strand")])
colnames(can) <- c("chromosome",
                   "external_gene_name",
                   "ensembl_gene_id",
                   "ensembl_transcript_id",
                   "transcription_start_site",
                   "strand")

ann_sub <- ann[category %in% c("Promoter",
                               "Enhancer - distal",
                               "Enhancer - proximal")][,c("chromosome",
                                                    "start",
                                                    "end",
                                                    "category",
                                                    "category_id",
                                                    "regulatory_feature_stable_id")]
colnames(ann_sub)[names(ann_sub) == "start"] <- "category_start"
colnames(ann_sub)[names(ann_sub) == "end"] <- "category_end"

### RUN

# define domains around tss
can_domain <- dominion(tss = can[,c("ensembl_transcript_id",
                                    "transcription_start_site",
                                    "strand")],
                       distal_extension = 300000)
can <- can[can_domain, on = c("ensembl_transcript_id",
                              "transcription_start_site",
                              "strand")]

# calculate tss domainoverlap for annotations 
      can_ann <- all_canonical(reg = ann_sub[,c("chromosome",
                                                "category_start",
                                                "category_end",
                                                "category_id",
                                                "regulatory_feature_stable_id")],
                               tss = can[, c("chromosome",
                                             "ensembl_transcript_id",
                                             "transcription_start_site",
                                             "domain",
                                             "domain_start",
                                             "domain_end" )],
                               chr = can$chromosome[1])
      can_ann <- ann_sub[can_ann, on = c("chromosome",
                                         "category_start",
                                         "category_end",
                                         "category_id",
                                         "regulatory_feature_stable_id")]
      out_all <- can_ann

# priorotise closest feature
out_1 <- ddply(out_all, "category", closest_one)

return(list(out_1, out_all))
}

### ARGUMENTS 

m.can.file = "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_canPC_pos.csv"
m.ann.file = "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
h.can.file = "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_canPC_pos.csv"
h.ann.file = "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/h.sap_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"

### IMPORT

m_can <- fread(m.can.file)
m_ann <- fread(paste0("gunzip -cq ", m.ann.file))
h_can <- fread(h.can.file)
h_ann <- fread(paste0("gunzip -cq ", h.ann.file))

### RUN

m_sod1 <- onix(can = m_can, ann = m_ann, gene_name = "Sod1")
h_sod1 <- onix(can = h_can, ann = h_ann, gene_name = "SOD1")

m_sod1_utr <- m_ann[transcript_id == m_can$ensembl_transcript_id[m_can$external_gene_name == "Sod1"][1] &
                      category == "Exon - UTR"][, 1:6]
h_sod1_utr <- h_ann[transcript_id == h_can$ensembl_transcript_id[h_can$external_gene_name == "SOD1"][1] &
                      category == "Exon - UTR"][, 1:6]
colnames(m_sod1_utr) <- c("chromosome", "category_start", "category_end","category", "category_id", "regulatory_feature_stable_id")
colnames(h_sod1_utr) <- c("chromosome", "category_start", "category_end","category", "category_id", "regulatory_feature_stable_id")

m_sod1_1 <- as.data.table(rbind.fill(m_sod1[[1]], m_sod1_utr))
h_sod1_1 <- as.data.table(rbind.fill(h_sod1[[1]], h_sod1_utr))



align <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/hsap_grch38_v_mmus_grcm38_v101_alignment.bed.gz")
h_align <- align[chromosome_human == h_sod1_1$chromosome[1]]
m_align <- align[chromosome_mouse == m_sod1_1$chromosome[1]]


# for each human annotation
# find alignment
# find conservation
# find coservation with orthologous element

out_list <- list()
for (j in 1:length(unique(h_sod1_1$category))){

h_sod1_sub2 <- h_sod1_1[category == unique(h_sod1_1$category)[j]]

# for each human annotation in chromosome:
n_total <- rep(NA, nrow(h_sod1_sub2))
n_aligned <- rep(NA, nrow(h_sod1_sub2))
n_conserved <- rep(NA, nrow(h_sod1_sub2))
n_conserved_2 <- rep(NA, nrow(h_sod1_sub2))

for(i in 1:nrow(h_sod1_sub2)){
  
  h_sod1_sub <- h_sod1_sub2[i,]
  
  n_total[i] <- sum(abs( (h_sod1_sub$category_end + 1) - h_sod1_sub$category_start))
  
  # orthologous sequences for human annotation to get total alignment
  ann_align <- orthologous.seq(alignment = h_align[,c("chromosome_human", "start_human", "end_human",
                                                    "chromosome_mouse", "start_mouse", "end_mouse")],
                               bed = h_sod1_sub[,c("chromosome", "category_start", "category_end")])
  n_aligned[i] <- sum(abs( (ann_align$end_A + 1) - ann_align$start_A))
  
  # total of mouse annotation in orthologous sequences for conservation
  m_ann_conserved <- bed.intersect(bed1 = ann_align[,c("chromosome_B", "start_B", "end_B")],
                                   bed2 = m_ann[category == h_sod1_sub$category[1]][,c("chromosome", "start", "end")])
  n_conserved[i] <- sum(abs ((m_ann_conserved$end + 1) - m_ann_conserved$start))
  
  m_ann_conserved_2 <- bed.intersect(bed1 = ann_align[,c("chromosome_B", "start_B", "end_B")],
                                   bed2 = m_sod1_1[category == h_sod1_sub$category[1]][,c("chromosome", "category_start", "category_end")])
  n_conserved_2[i] <- sum(abs ((m_ann_conserved_2$end + 1) - m_ann_conserved_2$start))
  
  print(i)
  }

out_list[[j]] <- data.table(category_id = h_sod1_sub2$category_id,
                          n_total = n_total,
                          n_aligned = n_aligned,
                          n_conserved = n_conserved,
                          n_conserved_2 = n_conserved_2)

}

h_out <- do.call("rbind", out_list)
h_output <- h_out[h_sod1_1, on = "category_id"]
h_output$ensembl_transcript_id <- h_output$ensembl_transcript_id[1]

m_output <- m_sod1_1
m_output$ensembl_transcript_id <- m_output$ensembl_transcript_id[1]

### EXPORT

fwrite(h_output, "~/Dropbox/h.sap.vs.m.mus_Ensembl_v101_SOD1_feature_alignment.csv")
fwrite(m_output, "~/Dropbox/m.mus_Ensembl_v101_Sod1_features.csv")


