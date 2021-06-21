### SCRIPT that qced and calculated constraint scores for windows

rm(list = ls())
graphics.off()

library(data.table)
library(MASS)

### SET ARGS 
contraint.variables.path <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/Variables/mmus_GRC38_constraint_var_by_window_1000_100_WMallSP_chr"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/mmus_GRC38_constraint_by_window_1000_100_WMallSP.csv.gz"
sm.max <- 500
cm.max <- 225
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

# # filter by n SNV percentile
# percentile <- ecdf(dt$n_SNV_weighted)
# n_SNV_percentile <- ceiling(percentile(dt$n_SNV_weighted)*100)
# rm.id <- c(which(n_SNV_percentile <= 5), which(n_SNV_percentile >= 95))
# if (length(rm.id) != 0){
#   removed <- rbind(removed, dt[rm.id,])
#   dt <- dt[-rm.id,]
# }

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

# fit lm
# mod <- lm(n_SNV_weighted ~ k7_psnv_weighted + n_CG_weighted + n_CG_unmeth_weighted + n_sm_weighted + n_cm_weighted, data = dt)
mod <- lm(n_SNV_weighted ~ n_CG_unmeth_weighted + n_sm_weighted + n_cm_weighted, data = dt)
summary(mod)

# model coeff
mod_summary <- as.data.table(summary(mod)$coef)
mod_summary$Covar <- row.names(summary(mod)$coef)
mod_summary$Adj_r2 <- summary(mod)$adj.r.squared

# residual score and percentile
dt$residual <- studres(mod)
percentile <- ecdf(dt$residual)
dt$residual_percentile <- ceiling(percentile(dt$residual)*100)

# percentile <- ecdf(dt$n_SNV_weighted)
# dt$residual_percentile <- ceiling(percentile(dt$n_SNV_weighted)*100)

out <- dt[, c("chromosome", "start", "end", "residual_percentile")]

### EXPORT
fwrite(out, out.file, compress = "gzip")


#####

# how does OE compare with residuals? 
# could non-linear model be used?




