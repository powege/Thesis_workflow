rm(list = ls())
graphics.off()

library(data.table)
library(MASS)

### IMPORT
dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_constraint_var_by_region_MGP_allMUSMUS_QCed.csv")

# fit lm 
mod <- lm(length_cm_weighted ~ f_CG, data = dt)
# summary(mod)
# plot(mod)

# model coeff
mod_summary <- as.data.table(summary(mod)$coef)
mod_summary$Covar <- row.names(summary(mod)$coef)
mod_summary$Adj_r2 <- summary(mod)$adj.r.squared
# error_95 <- qnorm(0.975)*sd(mod$residuals)/sqrt(length(mod$residuals))   # 95% confidence interval

# residual score
dt$residual <- studres(mod)

out <- dt[, c("chromosome", "start", "end", "residual")]

### EXPORT
fwrite(out, "~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_WG_constraint_by_region_MGP_allMUSMUS.csv")
fwrite(mod_summary, "~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_WG_constraint_by_region_MGP_allMUSMUS_coef.csv")




