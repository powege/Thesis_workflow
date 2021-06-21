rm(list=ls())
graphics.off()

library(data.table)
library(pROC)

### FUNCTIONS

# OMIM = vSTATUS[1]
# df <- dt_viability
kapow <- function(OMIM, VARS, df){
  
  x <- c(OMIM, VARS)
  data <- df[, x, with = F]
  
  variable <- names(data)[-1]
  group <- rep(OMIM, length(variable))
  Estimate.vec <- rep(NA, length(variable))
  Std_Err.vec <- rep(NA, length(variable))
  P.vec <- rep(NA, length(variable))
  ROC.vec <- rep(NA, length(variable))
  N.vec <- rep(NA, length(variable))
  
  for(i in seq_along(variable)){
    
    ind <- !is.na(data[,(i+1), with = F])[,1]
    glm.data <- data[ind,]
    mod <- glm(reformulate(variable[i], names(glm.data)[1]), data = glm.data, family=binomial)
    # summary(mod)
    
    N.vec[i] <- table(glm.data[,1])[2]
    Estimate.vec[i] <- summary(mod)$coefficients[2,1]
    Std_Err.vec[i] <- summary(mod)$coefficients[2,2]
    P.vec[i] <- summary(mod)$coefficients[2,4]
    
    pred <- predict(mod, glm.data)
    obs <- as.numeric(unlist(glm.data[,1]))
    ROC.vec[i] <- auc(roc(obs, pred))
  }
  
  out <- data.table(group = group,
                    variable = variable, 
                    n = N.vec,
                    estimate = Estimate.vec,
                    estimate_se = Std_Err.vec,
                    p_val = P.vec,
                    ROC = ROC.vec)
  
  return(out)
}

# OMIM = "L_SV"
# df = dt_viability
kapow2 <- function(OMIM, VARS, df){
  
  x <- c(OMIM, VARS)
  glm.data <- df[, x, with = F]
  col10 <- names(glm.data)[-1]
  
  Pred.list <- list()
  for(i in seq_along(col10)){
    mod <- glm(reformulate(col10[i], names(glm.data)[1]), data = glm.data, family=binomial)
    pred <- predict(mod, glm.data)
    obs <- as.numeric(unlist(glm.data[,1]))
    Pred.list[[i]] <- data.frame(mmus_external_gene_name = df$mmus_external_gene_name,
                                 obs = obs,
                                 pred = pred)
    colnames(Pred.list[[i]]) <- c("mmus_external_gene_name", "obs", paste(col10[i], colnames(Pred.list[[i]])[3], sep = "_"))
  }
  
  return(Pred.list)
}

### SET VARS
infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv"
ROC.out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_constraint_score_IMPC_viability_fertility_ROC.csv"
ROC.raw.out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_constraint_score_IMPC_viability_prediction.csv"

### IMPORT
dt <- fread(infile)

### FORMAT
dt_viability <- dt[,c("mmus_external_gene_name",
                      "Z_nonsynonymous",
                      "OE_nonsynonymous",
                      "rat_dNdS",
                      "spretus_dNdS",
                      "pLI",
                      "oe_lof_upper",
                      "mis_z",
                      # "Z_nonsynonymous_percentile",
                      # "OE_nonsynonymous_percentile",
                      # "rat_dNdS_percentile",
                      # "spretus_dNdS_percentile",
                      # "pLI_percentile",
                      # "oe_lof_upper_percentile",
                      # "mis_z_percentile",
                      "viability_status")]
# dt_viability <- dt_viability[complete.cases(dt_viability) & viability_status != "",]
dt_viability <- dt_viability[viability_status != "",]
dt_fertility <- dt[,c("mmus_external_gene_name",
                      "Z_nonsynonymous",
                      "OE_nonsynonymous",
                      "rat_dNdS",
                      "spretus_dNdS",
                      "pLI",
                      "oe_lof_upper",
                      "mis_z",
                      # "Z_nonsynonymous_percentile",
                      # "OE_nonsynonymous_percentile",
                      # "rat_dNdS_percentile",
                      # "spretus_dNdS_percentile",
                      # "pLI_percentile",
                      # "oe_lof_upper_percentile",
                      # "mis_z_percentile",
                      "fertility_status")]
# dt_fertility <- dt_fertility[complete.cases(dt_fertility) & fertility_status != "",]
dt_fertility <- dt_fertility[fertility_status != "",]

# percentiles
percentile <- ecdf(dt_viability$Z_nonsynonymous[!duplicated(dt_viability$mmus_external_gene_name)])
dt_viability$Z_nonsynonymous_percentile <- ceiling(percentile(dt_viability$Z_nonsynonymous)*100)
percentile <- ecdf(dt_viability$OE_nonsynonymous[!duplicated(dt_viability$mmus_external_gene_name)])
dt_viability$OE_nonsynonymous_percentile <- 101 - (ceiling(percentile(dt_viability$OE_nonsynonymous)*100))
percentile <- ecdf(dt_viability$spretus_dNdS[!duplicated(dt_viability$mmus_external_gene_name)])
dt_viability$spretus_dNdS_percentile <- 101 - (ceiling(percentile(dt_viability$spretus_dNdS)*100))
percentile <- ecdf(dt_viability$rat_dNdS[!duplicated(dt_viability$mmus_external_gene_name)])
dt_viability$rat_dNdS_percentile <- 101 - (ceiling(percentile(dt_viability$rat_dNdS)*100))
percentile <- ecdf(dt_viability$mis_z[!duplicated(dt_viability$mmus_external_gene_name)])
dt_viability$mis_z_percentile <- ceiling(percentile(dt_viability$mis_z)*100)
percentile <- ecdf(dt_viability$pLI[!duplicated(dt_viability$mmus_external_gene_name)])
dt_viability$pLI_percentile <- ceiling(percentile(dt_viability$pLI)*100)
percentile <- ecdf(dt_viability$oe_lof_upper[!duplicated(dt_viability$mmus_external_gene_name)])
dt_viability$oe_lof_upper_percentile <- 101 - (ceiling(percentile(dt_viability$oe_lof_upper)*100))

percentile <- ecdf(dt_fertility$Z_nonsynonymous[!duplicated(dt_fertility$mmus_external_gene_name)])
dt_fertility$Z_nonsynonymous_percentile <- ceiling(percentile(dt_fertility$Z_nonsynonymous)*100)
percentile <- ecdf(dt_fertility$OE_nonsynonymous[!duplicated(dt_fertility$mmus_external_gene_name)])
dt_fertility$OE_nonsynonymous_percentile <- 101 - (ceiling(percentile(dt_fertility$OE_nonsynonymous)*100))
percentile <- ecdf(dt_fertility$spretus_dNdS[!duplicated(dt_fertility$mmus_external_gene_name)])
dt_fertility$spretus_dNdS_percentile <- 101 - (ceiling(percentile(dt_fertility$spretus_dNdS)*100))
percentile <- ecdf(dt_fertility$rat_dNdS[!duplicated(dt_fertility$mmus_external_gene_name)])
dt_fertility$rat_dNdS_percentile <- 101 - (ceiling(percentile(dt_fertility$rat_dNdS)*100))
percentile <- ecdf(dt_fertility$pLI[!duplicated(dt_fertility$mmus_external_gene_name)])
dt_fertility$mis_z_percentile <- ceiling(percentile(dt_fertility$mis_z)*100)
percentile <- ecdf(dt_fertility$pLI[!duplicated(dt_fertility$mmus_external_gene_name)])
dt_fertility$pLI_percentile <- ceiling(percentile(dt_fertility$pLI)*100)
percentile <- ecdf(dt_fertility$oe_lof_upper[!duplicated(dt_fertility$mmus_external_gene_name)])
dt_fertility$oe_lof_upper_percentile <- 101 - (ceiling(percentile(dt_fertility$oe_lof_upper)*100))

# set bernouli variables
dt_viability$L_SV <- 0
dt_viability$L_SV[dt_viability$viability_status %in% c("L", "SV")] <- 1
dt_viability$VN <- 0
dt_viability$VN[dt_viability$viability_status %in% c("VN")] <- 1
dt_fertility$Infertility <- 0
dt_fertility$Infertility[dt_fertility$fertility_status %in% c("FI", "FIMI", "MI")] <- 1

## calculate roc

# remove NAs
vSTATUS = c(
  "L_SV",
  "VN"
)
fSTATUS = c(
  "Infertility"
)
VARS <- c(
  "Z_nonsynonymous_percentile",
  "OE_nonsynonymous_percentile",
  # "rat_dNdS_percentile",
  "spretus_dNdS_percentile",
  "pLI_percentile",
  "oe_lof_upper_percentile",
  "mis_z_percentile"
)
v.out.list <- list()
for (i in 1:length(vSTATUS)){
  v.out.list[[i]] <- kapow(vSTATUS[i], VARS, dt_viability)
}
f.out.list <- list()
for (i in 1:length(fSTATUS)){
  f.out.list[[i]] <- kapow(fSTATUS[i], VARS, dt_fertility)
}
ROC_table <- rbind( do.call("rbind", v.out.list), do.call("rbind", f.out.list) )

# get prediction for lethality
predict_raw <- kapow2("L_SV", VARS, dt_viability)
predict_raw <- do.call("cbind", predict_raw)
predict_raw <- predict_raw[, !duplicated(colnames(predict_raw))]


### EXPORT 
fwrite(ROC_table, ROC.out.file)
fwrite(predict_raw, ROC.raw.out.file)

#####

out2 <- ROC_table[variable %in% c("OE_nonsynonymous_percentile", "spretus_dNdS_percentile") & group %in% c("L_SV")]
out2$estimate <- round(out2$estimate, 2)
out2$estimate_se <- round(out2$estimate_se, 2)
out2$p_val <- round(out2$p_val, 2)
out2$ROC <- round(out2$ROC, 2)
fwrite(out2, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Figures_tables/table_supplementary_constraint_score_IMPC_viability_ROC.csv")


####

# library(data.table)
# library(ggplot2)
# library(plotROC)
# library(pROC)
# library("cowplot")
# 
# colnames(ROC_df1) <- c("H_external_gene_name", "obs", 
#                        paste0("missense Z-score\nAUC=", auc.mizZ, " p=", p.mizZ),
#                        paste0("pLI\nAUC=", auc.pLI, " p=", p.pLI),                  
#                        paste0("RVIS\nAUC=", auc.RVIS, " p=", p.RVIS),    
#                        paste0("Human funZ\nAUC=", auc.HfunZ, " p=", p.HfunZ),
#                        paste0("Mouse funZ\nAUC=", auc.MfunZ, " p=", p.MfunZ))
# ROC_df1 <- melt(ROC_df1, id.vars = c("H_external_gene_name", "obs"))
# ROC_df1$variable <- as.factor(ROC_df1$variable)
# unique(ROC_df1$variable)
# ROC_df1$variable <- factor(ROC_df1$variable, levels = c(
#   "pLI\nAUC=0.78 p=9.66e-31", 
#   "Mouse funZ\nAUC=0.77 p=2.21e-25",
#   "Human funZ\nAUC=0.77 p=4.52e-26",
#   "missense Z-score\nAUC=0.75 p=9.17e-23", 
#   "RVIS\nAUC=0.74 p=3.34e-22"
# ))
# 
# plotB <- ggplot(ROC_df, aes(d = obs, m = value, color = variable)) + 
#   geom_roc(n.cuts = 0) + 
#   style_roc() +
#   theme(legend.key.size = unit(1, "cm"),
#         legend.title = element_blank(),
#         legend.position = c(0.75, 0.25),
#         text = element_text(size = 14),
#         panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.background=element_blank(),
#         plot.margin=unit(c(1,1,1,1),"cm"))
# plotB

