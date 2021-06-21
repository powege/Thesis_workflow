rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS

### FUNCTION that returns the seqence match between bed file and vector
# (ie for each sequence in bed file1, how many integers overlap with vector?)
bed_match <- function(bed, vec){
  
  require(data.table)
  
  bed_dt <- setDT(bed) # set as bed as data.table
  colnames(bed_dt) <- c("start", "end")
  bed_dt[, ind := .I] # add uniqe index to data.table
  setkey(bed_dt) # sets keys // order data by all columns
  
  vec_dt <- as.data.table(vec, key = 'vec') # convert to data.table
  vec_dt[, vec2 := vec] # dublicate column
  
  # Fast overlap join:
  ans1 = foverlaps(vec_dt, bed_dt, by.x = c('vec', 'vec2'), by.y = c('start', 'end'),
                   type = "within", nomatch = 0L)
  counts <- ans1[, .N, keyby = ind] # count by ind
  # merge to inital data
  bed_dt[, n_match := counts[bed_dt, on = .(ind), x.N]]
  setorder(bed_dt, ind) # reorder by ind to get inital order
  bed_dt[, ind := NULL] # deletes ind colum
  bed_dt[is.na(n_match), n_match := 0L] # NAs is 0 count
  
  return(bed_dt$n_match) 
}

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

#### SET PATHS
pos.in.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_canPC_pos.csv"
meth.in.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/m.mus_grc38_ensembl_transcript_methylation.csv"

### IMPORT 
pos <- fread(pos.in.file)
meth_all <- fread(paste0("gunzip -cq ", meth.in.file))

### FORMAT

pos <- pos[!is.na(genomic_coding_start) & !is.na(genomic_coding_end),] # remove exons with no CDS POS 

# number of CGs
pos_list <- list()
for (chr in 1:19){
  pos_sub <- pos[chromosome_name == chr]
  pos_sub$n_CG <- bed_match(bed = pos_sub[,c("genomic_coding_start", "genomic_coding_end")],
                                   vec = c(meth_all$start[meth_all$chromosome == chr], meth_all$end[meth_all$chromosome == chr]))
  pos_list[[chr]] <- pos_sub
  print(chr)
}
pos <- do.call("rbind", pos_list)

# number of unmethylated CGs
unmeth <- meth_all[coverage >= 5 & percent_methylated <= 20][,c("chromosome", "start", "end")]
pos_list <- list()
for (chr in 1:19){
  pos_sub <- pos[chromosome_name == chr]
  pos_sub$n_CGunmethylated <- bed_match(bed = pos_sub[,c("genomic_coding_start", "genomic_coding_end")],
                            vec = c(unmeth$start[unmeth$chromosome == chr], unmeth$end[unmeth$chromosome == chr]))
  pos_list[[chr]] <- pos_sub
  print(chr)
}
pos <- do.call("rbind", pos_list)

# number of methylated CGs
meth <- meth_all[coverage > 5 & percent_methylated > 60][,c("chromosome", "start", "end")]
pos_list <- list()
for (chr in 1:19){
  pos_sub <- pos[chromosome_name == chr]
  pos_sub$n_CGmethylated <- bed_match(bed = pos_sub[,c("genomic_coding_start", "genomic_coding_end")],
                                        vec = c(meth$start[meth$chromosome == chr], meth$end[meth$chromosome == chr]))
  pos_list[[chr]] <- pos_sub
  print(chr)
}
pos <- do.call("rbind", pos_list)

# aggregate by transcript
output <- aggregate(list(n_CG = pos$n_CG, n_CGmethylated = pos$n_CGmethylated, n_CGunmethylated = pos$n_CGunmethylated), 
                    by=list(ensembl_transcript_id=pos$ensembl_transcript_id), FUN=sum)

### EXPORT 
fwrite(output, out.file)





