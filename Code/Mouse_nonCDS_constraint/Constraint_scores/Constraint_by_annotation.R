### SCRIPT that qced and calculated constraint scores by annotation

rm(list = ls())
graphics.off()

library(data.table)
library(MASS)

### SET ARGS 
constraint.variables.file <- "/Documents and Settings/ASUS/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/Variables/mmus_GRC38_constraint_var_by_annotation_WMallSP.csv.gz"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/mmus_GRC38_constraint_by_annotation_WMallSP.csv.gz"
sm.max <- 1
cm.max <- 0.1

### IMPORT

dt <- fread(paste0("gunzip -cq ", constraint.variables.file))

### FORMAT

dt$length <- ( dt$end +1 ) - dt$start

### QC

removed <- data.frame() # set dt for removed 

# filter by N mask fraction
rm.id <- which(dt$n_Nm > 0)
if (length(rm.id) != 0){
  removed <- rbind(removed, dt[rm.id,])
  dt <- dt[-rm.id,]
}

# filter by soft mask fraction
rm.id <- which(dt$n_sm >= sm.max)
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




