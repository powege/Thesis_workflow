rm(list = ls())
graphics.off()

library(data.table)

# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) { stop("arguments must be supplied", call.=FALSE) } 

# Correlation between regulatory feature constraint (UTRs, promoters, enhancers, CTCF binding sites) 
# and GERP conservation.

### FUNCTIONS

# FUNCTION that returns the seqence match between bed file and vector
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

### FUNCTION vectorisation of seq
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### FUNCTION
# cs <- cs
# score <- "OER_percentile"
# ann_sub <- ann[annotation == "Exon - UTR"]
# n_chr <- 19
# slide_size <- 50
alakazam <- function(cs, score, gerp, n_chr){
  
  # subset score of interest
  cs_sub <- cs[,c("chromosome", "start", "end", "annotation", score), with=F]
  colnames(cs_sub) <- c("chromosome", "start", "end", "annotation", "percentile")
  
  chr_list <- list()
  for (chr in 1:n_chr){
    
    cs_chr <- cs_sub[chromosome == chr]
    bed <- cs_chr[,c("start", "end")]
    
    ann_chr <- gerp[chromosome == chr]
    ann_vec <- as.vector(unique(unlist(seq2(from = ann_chr$start, ann_chr$end))))
    cs_chr$match <- bed_match(bed, ann_vec)
    
    chr_list[[chr]] <- cs_chr
    print(chr)
  }
  dt_sub <- do.call("rbind", chr_list)
  dt_sub$fraction <- dt_sub$match / ((dt_sub$end + 1) - dt_sub$start)
  return(dt_sub)
}


### RUN FOR MOUSE

### SET ARGS

# cs_file <- args[1]
# ann_file <- args[2]
# score <- as.character(args[3])
# n_chr <- as.numeric(args[4])
# out_file <- args[5]

cs_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_CS.csv"
gerp_file <- "~/Dropbox/PhD/Data/Ensembl/GERP/gerp_constrained_elements.mus_musculus.bed"
score <- "OE_ratio_percentile"
n_chr <- 19
out_file <- "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_GERP_by_annotation_percentile.csv"

### IMPORT

cs <- fread(cs_file)
gerp <- fread(gerp_file)

### FORMAT

# colnames
colnames(gerp) <- c("chromosome", "start", "end")

# percentiles
tmp_list <- list()
for (cat in 1:length(unique(cs$annotation))){
  tmp_sub <- cs[annotation == unique(cs$annotation)[cat]]
  
  percentile <- ecdf(tmp_sub$residual[!duplicated(tmp_sub$ID)])
  tmp_sub$residual_percentile <- percentile(tmp_sub$residual)
  
  percentile <- ecdf(tmp_sub$OE_ratio[!duplicated(tmp_sub$ID)])
  tmp_sub$OE_ratio_percentile <- percentile(tmp_sub$OE_ratio)
  
  tmp_list[[cat]] <- tmp_sub
}
cs <- do.call("rbind", tmp_list)
rm(tmp_list, tmp_sub)

# bin annotation constraint percentiles
cs$residual_percentile <- ceiling(cs$residual_percentile * 100)
cs$OE_ratio_percentile <- ceiling(cs$OE_ratio_percentile * 100)


### RUN

output <- alakazam(cs=cs, 
                   score = score, 
                   gerp = gerp, 
                   n_chr = n_chr)

### EXPORT
  
fwrite(output, out_file)



#####

# ### STACK OVERFLOW
# 
# library(data.table)
# set.seed(1)
# 
# dt1 <- data.table(start = sample(1:100000, 1000),
#                   end = sample(1:100000, 1000))
# dt2 <- data.table(start = sample(1:100000, 10),
#                   end = sample(1:100000, 10))
# dt1[ start > end, `:=`( start = end, end = start)] # ensure start <= end
# dt2[ start > end, `:=`( start = end, end = start)]
# 
# 
# seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
# vec <- unique(unlist(seq2(from = dt2$start, to = dt2$end)))
# 
# dt1[, ind := .I] # add uniqe index to data.table
# setkey(dt1) # sets keys // order data by all columns
# vec_dt <- as.data.table(vec, key = 'vec') # convert to data.table
# vec_dt[, vec2 := vec] # dublicate column
# # Fast overlap join:
# ans1 = foverlaps(vec_dt, dt1, by.x = c('vec', 'vec2'), by.y = c('start', 'end'),
#                  type = "within", nomatch = 0L)
# counts <- ans1[, .N, keyby = ind] # count by ind
# # merge to inital data
# dt1[, n_match := counts[dt1, on = .(ind), x.N]]
# setorder(dt1, ind) # reorder by ind to get inital order
# dt1[, ind := NULL] # deletes ind colum
# dt1[is.na(n_match), n_match := 0L] # NAs is 0 count








