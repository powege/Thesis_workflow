rm(list = ls())
graphics.off()

library(data.table)

# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) { stop("arguments must be supplied", call.=FALSE) } 

### SET ARGS

# cs_file <- args[1]
# ann_file <- args[2]
# score <- as.character(args[3])
# n_chr <- as.numeric(args[4])
# slide_size <- as.numeric(args[5])
# out_file <- args[6]

cs_file <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_constraint_by_window_750_50_MGP_allMUSMUS.csv"
snv_file <- "~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_mouse_mapping.csv"
score <- "OE"
n_chr <- 19
slide_size <- 50
out_file <- "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/ClinVar_SNV_by_OE_percentile_750_50_MGP_allMUSMUS.csv"

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
# score <- "OE_percentile"
# ann_sub <- ann[annotation == "Exon - UTR"]
# n_chr <- 19
# slide_size <- 50
alakazam <- function(cs, score, ann_sub, n_chr, slide_size){
  
  # subset score of interest
  cs_sub <- cs[,c("chromosome", "start", "end", score), with=F]
  colnames(cs_sub) <- c("chromosome", "start", "end", "percentile")
  
  chr_list <- list()
  for (chr in 1:n_chr){
    
    cs_chr <- cs_sub[chromosome == chr]
    bed <- cs_chr[,c("start", "end")]
    
    ann_chr <- ann_sub[chromosome == chr]
    ann_vec <- as.vector(unique(ann_chr$POS))
    cs_chr$match <- bed_match(bed, ann_vec)
    
    n_percentile <- as.data.table(table(cs_chr$percentile))
    colnames(n_percentile) <- c("percentile", "n_total")
    n_percentile$percentile <- as.numeric(n_percentile$percentile)
    n_percentile$n_total <- n_percentile$n_total * slide_size
    n_match <- aggregate(list(n_match=cs_chr$match), by=list(percentile=cs_chr$percentile), FUN=sum)
    output <- merge(n_percentile, n_match)
    output$chromosome <- chr
    
    chr_list[[chr]] <- output
    print(chr)
  }
  dt_sub <- do.call("rbind", chr_list)
  n_total <- aggregate(list(n_total=dt_sub$n_total), by=list(percentile=dt_sub$percentile), FUN=sum)
  n_match <- aggregate(list(n_match=dt_sub$n_match), by=list(percentile=dt_sub$percentile), FUN=sum)
  dt_sub <- merge(n_total, n_match)
  dt_sub$fraction <- dt_sub$n_match / dt_sub$n_total
  
  dt_sub <- dt_sub[,c("percentile", "fraction")]
  return(dt_sub)
}


### IMPORT

cs <- fread(cs_file)
snv <- fread(snv_file)

### FORMAT

snv <- snv[,c(6,7,10)]
colnames(snv) <- c("chromosome", "POS", "annotation")

# split coding and non-coding 
cds <- snv[annotation %like% "A"]
noncds <- snv[!annotation %like% "A"]

# percentiles
percentile <- ecdf(cs$residual)
cs$residual_percentile <-ceiling(percentile(cs$residual) * 100)
cs$residual <- NULL
percentile <- ecdf(cs$OE)
cs$OE_percentile <-ceiling(percentile(cs$OE) * 100)
cs$OE <- NULL
# percentile <- ecdf(cs$OER)
# cs$OER_percentile <-ceiling(percentile(cs$OER) * 100)
# cs$OER <- NULL

### RUN

cds_out <- alakazam(cs=cs, 
                    score = score, 
                    ann_sub = cds, 
                    n_chr = n_chr, 
                    slide_size = slide_size)
cds_out$category <- "CDS"
noncds_out <- alakazam(cs=cs, 
                    score = score, 
                    ann_sub = noncds, 
                    n_chr = n_chr, 
                    slide_size = slide_size)
noncds_out$category <- "nonCDS"

output <- rbind(cds_out, noncds_out)

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








