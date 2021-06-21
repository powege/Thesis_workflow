rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {  stop("At least one argument must be supplied", call.=FALSE) }


### SET ARGS 
contraint_variables_file <- args[1]
out_file <- args[2]
species <- args[3]
length_min <- as.numeric(args[4])
length_max <- as.numeric(args[5])
sm_max <- as.numeric(args[6])
cm_max <- as.numeric(args[7])

# contraint_variables_file <- "~/Dropbox/PhD/Data/NC_constraint/Variables/Human_constraint_var_by_region_gnomad_MAF001_chr1.csv"
# out_file <- ""
# species <- "human"
# length_min <- 2
# length_max <- 10000
# sm_max <- 1
# cm_max <- 0.5

# contraint_variables_file <- "/well/lindgren/George/Data/NC_constraint/Human_constraint_var_by_region_gnomad_MAF001_chr1.csv"
# out_file <- "/well/lindgren/George/Data/NC_constraint/Human_constraint_var_by_region_gnomad_MAF001_QCed.csv"
# species <- "human"
# length_min <- 2
# length_max <- 10000
# sm_max <- 1
# cm_max <- 0.5

# if (species == "mouse") { chr_vec <- c(1:19)}
# if (species == "human") { chr_vec <- c(1:22)}


### IMPORT

dt <- fread(contraint_variables_file)

### FORMAT 

# colnames(dt) <- c("chromosome", "start", "end", "f_CG", "f_sm", "f_cm", "length_cm_weighted")
colnames(dt) <- c("chromosome","start","end","n_CG","f_CG","n_sm","f_sm","n_cm","f_cm","length","length_cm_weighted")

dt <- dt[,c("chromosome", "start", "end", "f_CG", "f_sm", "f_cm", "length_cm_weighted")]


### QC

removed <- data.frame() # set dt for removed 

# filter by complete cases
removed <- dt[!complete.cases(dt),]
dt <- dt[complete.cases(dt),]

# filter by length
# hist(dt$length)
rm.id <- which(dt$length < length_min | dt$length > length_max)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by soft mask fraction
# hist(dt$f_sm)
rm.id <- which(dt$f_sm >= sm_max)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by coverage mask fraction
# hist(dt$f_cm)
rm.id <- which(dt$f_cm >= (cm_max))
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

output <- dt[,c("chromosome", "start", "end", "f_CG", 'f_sm', "length_cm_weighted")]


### OUTPUT
fwrite(output, out_file, append = T, col.names = T)

#####

# x <- as.data.table(table(dt$length_cm_weighted))
# sub <- dt[which(dt$length_cm_weighted %% 1 != 0)]
# summary(lm(length_cm_weighted ~ f_CG + f_sm, data = dt))



