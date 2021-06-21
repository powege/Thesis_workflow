### SCRIPT that creates null datasets with same k11_from_N as observed:
# 1. where k11_to_N is sampled from a binomial distribution with same pMU as k1CG model
# ... k3 - k9
# 2. where k11_to_N is sampled from a binomial distribution with same pMU as k11 model 

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 

### SET VARS
k11.counts.chr.infile <- args[1]

k11.p.snv.infile <- args[2]
k9.p.snv.infile <- args[3]
k7.p.snv.infile <- args[4]
k5.p.snv.infile <- args[5]
k3.p.snv.infile <- args[6]
k1CG.p.snv.infile <- args[7]

k11.null.outfile <- args[8]
k9.null.outfile <- args[9]
k7.null.outfile <- args[10]
k5.null.outfile <- args[11]
k3.null.outfile <- args[12]
k1CG.null.outfile <- args[13]


### IMPORT 
k11_chr <- fread(paste0("gunzip -cq ", k11.counts.chr.infile))

k11_p_snv <- fread(paste0("gunzip -cq ", k11.p.snv.infile))
k9_p_snv <- fread(paste0("gunzip -cq ", k9.p.snv.infile))
k7_p_snv <- fread(paste0("gunzip -cq ", k7.p.snv.infile))
k5_p_snv <- fread(paste0("gunzip -cq ", k5.p.snv.infile))
k3_p_snv <- fread(paste0("gunzip -cq ", k3.p.snv.infile))
k1CG_p_snv <- fread(paste0("gunzip -cq ", k1CG.p.snv.infile))


### FORMAT

names(k11_chr)[names(k11_chr) == 'k11_to'] <- 'to'
k11_chr <- k11_chr[,c("k11_from", "to", "k11_from_N", "chromosome")]

k11_p_snv$k9_from <- stri_sub(k11_p_snv$k11_from, 2, -2)
k11_p_snv$k7_from <- stri_sub(k11_p_snv$k11_from, 3, -3)
k11_p_snv$k5_from <- stri_sub(k11_p_snv$k11_from, 4, -4)
k11_p_snv$k3_from <- stri_sub(k11_p_snv$k11_from, 5, -5)
k11_p_snv <- k9_p_snv[k11_p_snv, on = c("k9_from", "to")]
k11_p_snv <- k7_p_snv[k11_p_snv, on = c("k7_from", "to")]
k11_p_snv <- k5_p_snv[k11_p_snv, on = c("k5_from", "to")]
k11_p_snv <- k3_p_snv[k11_p_snv, on = c("k3_from", "to")]
k11_p_snv <- k1CG_p_snv[k11_p_snv, on = c("k3_from", "to")]

k11 <- unique(k11_p_snv[,c("k11_from", "to", "k11_mu_rate")])
k9 <- unique(k11_p_snv[,c("k11_from", "to", "k9_mu_rate")])
k7 <- unique(k11_p_snv[,c("k11_from", "to", "k7_mu_rate")])
k5 <- unique(k11_p_snv[,c("k11_from", "to", "k5_mu_rate")])
k3 <- unique(k11_p_snv[,c("k11_from", "to", "k3_mu_rate")])
k1CG <- unique(k11_p_snv[,c("k11_from", "to", "k1CG_mu_rate")])
rm(k11_p_snv)

k11$k11_mu_rate[k11$k11_mu_rate == 0] <- min(k11$k11_mu_rate[k11$k11_mu_rate > 0])
k9$k9_mu_rate[k11$k9_mu_rate == 0] <- min(k11$k11_mu_rate[k11$k11_mu_rate > 0])
k7$k7_mu_rate[k11$k7_mu_rate == 0] <- min(k11$k11_mu_rate[k11$k11_mu_rate > 0])
k5$k5_mu_rate[k11$k5_mu_rate == 0] <- min(k11$k11_mu_rate[k11$k11_mu_rate > 0])
k3$k3_mu_rate[k11$k3_mu_rate == 0] <- min(k11$k11_mu_rate[k11$k11_mu_rate > 0])

k11_null <- k11_chr[k11, on = c("k11_from", "to")]
k9_null <- k11_chr[k9, on = c("k11_from", "to")]
k7_null <- k11_chr[k7, on = c("k11_from", "to")]
k5_null <- k11_chr[k5, on = c("k11_from", "to")]
k3_null <- k11_chr[k3, on = c("k11_from", "to")]
k1CG_null <- k11_chr[k1CG, on = c("k11_from", "to")]

k11_null$k11_to_N <- rbinom(size = k11_null$k11_from_N, n = nrow(k11_null), prob = k11_null$k11_mu_rate)
k9_null$k11_to_N <- rbinom(size = k9_null$k11_from_N, n = nrow(k9_null), prob = k9_null$k9_mu_rate)
k7_null$k11_to_N <- rbinom(size = k7_null$k11_from_N, n = nrow(k7_null), prob = k7_null$k7_mu_rate)
k5_null$k11_to_N <- rbinom(size = k5_null$k11_from_N, n = nrow(k5_null), prob = k5_null$k5_mu_rate)
k3_null$k11_to_N <- rbinom(size = k3_null$k11_from_N, n = nrow(k3_null), prob = k3_null$k3_mu_rate)
k1CG_null$k11_to_N <- rbinom(size = k1CG_null$k11_from_N, n = nrow(k1CG_null), prob = k1CG_null$k1CG_mu_rate)

k11_null <- k11_null[,c("k11_from", "k11_from_N", "to", "k11_to_N", "chromosome")]
k9_null <- k9_null[,c("k11_from", "k11_from_N", "to", "k11_to_N", "chromosome")]
k7_null <- k7_null[,c("k11_from", "k11_from_N", "to", "k11_to_N", "chromosome")]
k5_null <- k5_null[,c("k11_from", "k11_from_N", "to", "k11_to_N", "chromosome")]
k3_null <- k3_null[,c("k11_from", "k11_from_N", "to", "k11_to_N", "chromosome")]
k1CG_null <- k1CG_null[,c("k11_from", "k11_from_N", "to", "k11_to_N", "chromosome")]

### EXPORT 
fwrite(k11_null, k11.null.outfile, compress = "gzip")
fwrite(k9_null, k9.null.outfile, compress = "gzip")
fwrite(k7_null, k7.null.outfile, compress = "gzip")
fwrite(k5_null, k5.null.outfile, compress = "gzip")
fwrite(k3_null, k3.null.outfile, compress = "gzip")
fwrite(k1CG_null, k1CG.null.outfile, compress = "gzip")






