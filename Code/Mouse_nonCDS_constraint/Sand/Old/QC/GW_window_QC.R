rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {  stop("At least one argument must be supplied", call.=FALSE) }


### SET ARGS 
contraint_variables_path <- args[1]
out_file <- args[2]
species <- args[3]
sm_max <- as.numeric(args[4])
cm_max <- as.numeric(args[5])

# contraint_variables_path <- "/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_MGP_allSTRAIN_chr"
# out_file <- "/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_MGP_allSTRAIN_QCed.csv"
# species <- "mouse"
# sm_max <- 0.9
# cm_max <-0.5

# contraint_variables_path <- "/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_wild_chr"
# out_file <- "/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_wild_QCed.csv"
# species <- "mouse"
# sm_max <- 0.9
# cm_max <-0.5

if (species == "mouse") { chr_vec <- c(1:19)}
if (species == "human") { chr_vec <- c(1:22)}


### IMPORT

# all chr
dt.list <- list()
for (chr in chr_vec){
  dt.list[[chr]] <- fread(paste0(contraint_variables_path, chr, ".csv"))
}
dt <- rbind.fill(dt.list)
rm(dt.list)

### QC

removed <- data.frame() # set dt for removed 

# filter by complete cases
removed <- dt[!complete.cases(dt),]
dt <- dt[complete.cases(dt),]

# filter by soft mask fraction
# hist(dt$n_sm)
# rm.id <- which(dt$n_sm >= (sm_max*750))
# if (length(rm.id) != 0){
#   removed <- rbind(removed, dt[rm.id,])
#   dt <- dt[-rm.id,]
# }
# hist(dt$n_sm_weighted)
rm.id <- which(dt$n_sm_weighted >= (sm_max*400))
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by coverage mask fraction
# hist(dt$n_cm)
# rm.id <- which(dt$n_cm >= (cm_max*750))
# if (length(rm.id) != 0){
#   removed <- rbind(removed, dt[rm.id,])
#   dt <- dt[-rm.id,]
# }
# hist(dt$n_cm_weighted)
rm.id <- which(dt$n_cm_weighted >= (cm_max*400))
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

### OUTPUT
fwrite(dt, out_file)

#####