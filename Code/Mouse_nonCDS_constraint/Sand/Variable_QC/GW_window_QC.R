rm(list = ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {  stop("At least one argument must be supplied", call.=FALSE) }

### SET ARGS 
contraint.variables.path <- args[1]
out.file <- args[2]
rm.out.file <- args[3]
sm.max <- as.numeric(args[4])
cm.max <- as.numeric(args[5])

contraint.variables.path <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/Variables/mmus_GRC38_constraint_var_by_winndow_750_50_WMallSP_chr"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/Variables/mmus_GRC38_constraint_var_by_window_750_50_WMallSP_QCed_chr19.csv.gz"
rm.out.file <- ""
sm.max <- 0.9
cm.max <- 0.5


### IMPORT

# all chr
dt_list <- list()
for (chr in 1:19){
  dt_list[[chr]] <- fread(paste0("gunzip -cq ", contraint.variables.path, chr, ".csv.gz"))
}
dt <- do.call("rbind", dt_list)
rm(dt_list)

## QC

## variable distributions
hist(dt$n_SNV_weighted)
hist(dt$k7_psnv_weighted)
hist(dt$n_sm_weighted)
hist(dt$n_Nm_weighted)
hist(dt$n_cm_weighted)

removed <- data.frame() # set dt for removed 

# filter by N mask fraction
rm.id <- which(dt$n_Nm_weighted > 0)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by soft mask fraction
rm.id <- which(dt$n_sm_weighted >= (sm.max*400))
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by coverage mask fraction
rm.id <- which(dt$n_cm_weighted >= (cm.max*400))
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by complete cases
removed <- rbind(removed, dt[!complete.cases(dt),])
dt <- dt[complete.cases(dt),]


### OUTPUT
fwrite(dt, out.file, compress = "gzip")
fwrite(removed, rm.out.file, compress = "gzip")


#####