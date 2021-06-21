### SCRIPT that qced and calculated constraint scores for windows

rm(list = ls())
graphics.off()

library(data.table)
library(MASS)

### SET ARGS 
contraint.variables.path <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/Variables/mmus_GRC38_constraint_var_by_window_1000_100_WMallSP_chr"
out.file <- ""
sm.max <- 500
cm.max <- 500
snv.min <- 5

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
# hist(dt$n_SNV_weighted)
# hist(dt$k7_psnv_weighted)
# hist(dt$n_unmeth_weighted)
# hist(dt$n_sm_weighted)
# hist(dt$n_Nm_weighted)
# hist(dt$n_cm_weighted)

removed <- data.frame() # set dt for removed 

# filter by n SNV percentile
percentile <- ecdf(dt$n_SNV_weighted)
n_SNV_percentile <- ceiling(percentile(dt$n_SNV_weighted)*100)
rm.id <- c(which(n_SNV_percentile <= 5), which(n_SNV_percentile >= 95))
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# # filter windows by n SNVs
# rm.id <- which(dt$n_SNV_weighted < snv.min)
# if (length(rm.id) != 0){
#   removed <- rbind(removed, dt[rm.id,])
#   dt <- dt[-rm.id,]
# }

# filter by N mask fraction
rm.id <- which(dt$n_Nm_weighted > 0)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by soft mask fraction
rm.id <- which(dt$n_sm_weighted >= sm.max)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by coverage mask fraction
rm.id <- which(dt$n_cm_weighted >= cm.max)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by complete cases
removed <- rbind(removed, dt[!complete.cases(dt),])
dt <- dt[complete.cases(dt),]

# calculate n_SNV_percentile
percentile <- ecdf(dt$n_SNV_weighted)
dt$n_SNV_percentile <- ceiling(percentile(dt$n_SNV_weighted)*100)

dt_md <- aggregate(list(avg_n_SNV=dt$n_SNV_weighted,
                       avg_k7_psnv=dt$k7_psnv_weighted,
                       avg_n_CG=dt$n_CG_weighted,
                       avg_n_unmeth=dt$n_CG_unmeth_weighted,
                       avg_n_sm=dt$n_sm_weighted,
                       avg_n_cm=dt$n_cm_weighted), by=list(SNV_percentile=dt$n_SNV_percentile), FUN=median)

dt_mn <- aggregate(list( avg_n_SNV=dt$n_SNV_weighted,
                         avg_k7_psnv=dt$k7_psnv_weighted,
                         avg_n_CG=dt$n_CG_weighted,
                         avg_n_unmeth=dt$n_CG_unmeth_weighted,
                         avg_n_sm=dt$n_sm_weighted,
                         avg_n_cm=dt$n_cm_weighted), by=list(SNV_percentile=dt$n_SNV_percentile), FUN=mean)

plot(dt_md$SNV_percentile, dt_md$avg_n_SNV)
plot(dt_md$SNV_percentile, dt_md$avg_k7_psnv)
plot(dt_md$SNV_percentile, dt_md$avg_n_CG)
plot(dt_md$SNV_percentile, dt_md$avg_n_unmeth)
plot(dt_md$SNV_percentile, dt_md$avg_n_sm)
plot(dt_md$SNV_percentile, dt_md$avg_n_cm)

plot(dt_mn$SNV_percentile, dt_mn$avg_n_SNV)
plot(dt_mn$SNV_percentile, dt_mn$avg_k7_psnv)
plot(dt_mn$SNV_percentile, dt_mn$avg_n_CG)
plot(dt_mn$SNV_percentile, dt_mn$avg_n_unmeth)
plot(dt_mn$SNV_percentile, dt_mn$avg_n_sm)
plot(dt_mn$SNV_percentile, dt_mn$avg_n_cm)

dt_1_10 <- dt[n_SNV_percentile %in% 1:10]
hist(dt_1_10$n_SNV_weighted)
hist(dt_1_10$k7_psnv_weighted)
hist(dt_1_10$n_CG_weighted)
hist(dt_1_10$n_CG_unmeth_weighted)
hist(dt_1_10$n_sm_weighted)
hist(dt_1_10$n_cm_weighted)

dt_11_20 <- dt[n_SNV_percentile %in% 11:20]
hist(dt_11_20$n_SNV_weighted)
hist(dt_11_20$k7_psnv_weighted)
hist(dt_11_20$n_CG_weighted)
hist(dt_11_20$n_CG_unmeth_weighted)
hist(dt_11_20$n_sm_weighted)
hist(dt_11_20$n_cm_weighted)



