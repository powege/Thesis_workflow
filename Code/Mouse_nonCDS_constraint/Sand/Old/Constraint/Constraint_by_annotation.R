rm(list = ls())
graphics.off()

library(data.table)
library(MASS)

### FUNCTIONS 
alakazam <- function(dt){
  
ann_vec <- unique(dt$annotation)
dt_list <- list()
summary_list <- list()

for (i in 1:length(ann_vec)){
  
  # subset annotaiton
  dt_sub <- dt[annotation == ann_vec[i]]
  
  # fit lm 
  mod <- lm(n_SNV ~ length + n_CG + n_sm, data = dt_sub)
  # summary(mod)
  # plot(mod)
  
  # model coeff
  tmp <- as.data.table(summary(mod)$coef)
  tmp$Covar <- row.names(summary(mod)$coef)
  tmp$Annotation <- ann_vec[i]
  tmp$Adj_r2 <- summary(mod)$adj.r.squared
  summary_list[[i]] <- tmp
  rm(tmp)
  # error_95 <- qnorm(0.975)*sd(mod$residuals)/sqrt(length(mod$residuals))   # 95% confidence interval

  # predict
  dt_sub$exp_SNV <- predict(mod, data.table(length = dt_sub$length,
                                            n_CG = dt_sub$n_CG,
                                            n_sm = dt_sub$n_sm))
  
  # residual score
  dt_sub$residual <- studres(mod)
  
  # OE ratio
  dt_sub$OE_ratio <- dt_sub$n_SNV / dt_sub$exp_SNV
  
  # OE ratio rank
  percentile <- ecdf(dt_sub$OE_ratio[!duplicated(dt_sub$ID)])
  dt_sub$OE_ratio_rank <- percentile(dt_sub$OE_ratio)
  dt_sub$OE_ratio_rank <- ceiling(dt_sub$OE_ratio_rank * 100)
  dt_sub <- dt_sub[order(OE_ratio_rank, -exp_SNV),]
  dt_sub$OE_ratio_rank <- 1:nrow(dt_sub)
  percentile <- ecdf(dt_sub$OE_ratio_rank[!duplicated(dt_sub$ID)])
  dt_sub$OE_ratio_rank <- percentile(dt_sub$OE_ratio_rank)
  
  # what is the predicted length for a sequence with one SNV?
  # mod_len <- lm(length ~ n_SNV, data = dt_sub)
  # tmp <- predict(mod_len, data.table(n_SNV = 1))
  
  # add to list
  dt_list[[i]] <- dt_sub
}

dt_out <- do.call("rbind", dt_list)
summary_out <- do.call("rbind", summary_list)
output <- list(dt_out, summary_out)

return(output)
}

### IMPORT
M_dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Variables/Mouse_GRC38_GENCODE_RegBuild_annotation_varQCed.csv")
H_dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Variables/Human_GRC38_GENCODE_RegBuild_annotation_varQCed.csv")

### RUN
M_dt_list <- alakazam(M_dt)
H_dt_list <- alakazam(H_dt)

M_dt_out <- M_dt_list[[1]]
H_dt_out <- H_dt_list[[1]]

M_summary_out <- M_dt_list[[2]]
H_summary_out <- H_dt_list[[2]]

### EXPORT
fwrite(M_dt_out, "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_CS.csv")
fwrite(H_dt_out, "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_CS.csv")

fwrite(M_summary_out, "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_CSmod.csv")
fwrite(H_summary_out, "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_CSmod.csv")



