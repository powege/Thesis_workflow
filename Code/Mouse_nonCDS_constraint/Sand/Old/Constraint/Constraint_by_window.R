rm(list = ls())
graphics.off()

library(data.table)
library(MASS)

### IMPORT
dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_constraint_var_by_window_750_50_wild_QCed.csv")

#####
dt <- dt[n_cm_weighted < (0.1*400)]
dt <- dt[n_sm_weighted < (1*400)]
#####

# fit lm 
mod <- lm(n_SNV_weighted ~ n_CG_weighted + n_sm_weighted, data = dt)
summary(mod)
# plot(mod)

# model coeff
mod_summary <- as.data.table(summary(mod)$coef)
mod_summary$Covar <- row.names(summary(mod)$coef)
mod_summary$Adj_r2 <- summary(mod)$adj.r.squared
# error_95 <- qnorm(0.975)*sd(mod$residuals)/sqrt(length(mod$residuals))   # 95% confidence interval

# # predict
dt$exp_SNV <- predict(mod, data.table(n_CG_weighted = dt$n_CG_weighted,
                                      n_sm_weighted = dt$n_sm_weighted,
                                      n_cm_weighted = dt$n_cm_weighted))

# residual score
dt$residual <- studres(mod)

# OE ratio
dt$OE <- dt$n_SNV_weighted / dt$exp_SNV

# # OE ratio rank
# percentile <- ecdf(dt$OE)
# dt$OER <- percentile(dt$OE)
# dt$OER <- ceiling(dt$OER * 100)
# dt <- dt[order(OER, -exp_SNV),]
# dt$OER <- 1:nrow(dt)
# percentile <- ecdf(dt$OER)
# dt$OER <- percentile(dt$OER)

out <- dt[, c("chromosome", "start", "end", 
              "residual",
              # "OER", 
              "OE"
              )]

### EXPORT
fwrite(out, "~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_constraint_by_window_750_50_wild.csv")
fwrite(mod_summary, "~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_constraint_by_window_750_50_wild_coef.csv")


#####

# # percentiles
# percentile <- ecdf(dt$residual)
# dt$residual_percentile <-ceiling(percentile(dt$residual) * 100)
# dt$residual <- NULL
# 
# p1 <- dt[residual_percentile == 1]
# p2 <- dt[residual_percentile == 2]
# p1to5 <- dt[residual_percentile %in% 1:5]
# 
# p100 <- dt[residual_percentile == 100]
# 
# hist(p1to5$n_SNV_weighted)
# hist(p1to5$n_CG_weighted)
# hist(p1to5$n_cm_weighted)
# hist(p1to5$n_sm_weighted)
# 
# hist(p2$n_SNV_weighted)
# hist(p2$n_CG_weighted)
# hist(p2$n_cm_weighted)
# hist(p2$n_sm_weighted)


