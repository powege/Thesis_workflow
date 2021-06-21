rm(list = ls())
graphics.off()

library(data.table)

### SET VARS
counts.all.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_CG_SNV_counts_by_chr.csv"
counts.meth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_CG_methylated_SNV_counts_by_chr.csv"
counts.unmeth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_SNV_counts_by_chr.csv"
mae.meth.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_obs/WM_Harr_etal_2016_allSPECIES_CG_resampling_by_chr_model_fit_specific_CG_methylated.csv"
mae.unmeth.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_obs/WM_Harr_etal_2016_allSPECIES_CG_resampling_by_chr_model_fit_specific_CG_unmethylated.csv"


### IMPORT
counts_all <- fread(counts.all.infile)
counts_m <- fread(counts.meth.infile)
counts_um <- fread(counts.unmeth.infile)


zapdos <- function(counts_all, counts_m, n_resample){
  
out_list <- list()
for(i in 1:n_resample){

# split into train and test sets
chr <- sample(1:19, 4)
train_null <- counts_all[!chromosome %in% chr] 
train <- counts_m[!chromosome %in% chr]
test <- counts_m[chromosome %in% chr]
train_null <- setDT(train_null)[, .(from_N = sum(from_N), to_N = sum(to_N)), by = .(from, to)] # sum across chromosomes
train <- setDT(train)[, .(from_N = sum(from_N), to_N = sum(to_N)), by = .(from, to)] # sum across chromosomes
test <- setDT(test)[, .(from_N = sum(from_N), to_N = sum(to_N)), by = .(from, to)] # sum across chromosomes

# calculate pmu
train_null$mu_rate <- train_null$to_N / train_null$from_N
train$mu_rate <- train$to_N / train$from_N
test$mu_rate <- test$to_N / test$from_N

# calculate AE
test <- test[,c("from", "to", "mu_rate")]
train_null <- train_null[,c("from", "to", "mu_rate")]
train <- train[,c("from", "to", "mu_rate")]
names(test)[names(test) == 'mu_rate'] <- 'mu_rate_obs'
names(train_null)[names(train_null) == 'mu_rate'] <- 'mu_rate_null'
names(train)[names(train) == 'mu_rate'] <- 'mu_rate_meth'
dt <- test[train, on = c("from", "to")]
dt <- dt[train_null, on = c("from", "to")]
dt$AE_null <- abs(dt$mu_rate_null - dt$mu_rate_obs)
dt$AE <- abs(dt$mu_rate_meth - dt$mu_rate_obs)
dt$test_chromosome <- paste(chr, collapse = "-")
out_list[[i]] <- dt[,c("from", "to", "AE_null", "AE", "test_chromosome")]
}
output <- do.call("rbind", out_list)
output$mutation <- paste0(output$from, ">", output$to)
return(output)
}

dt_meth <- zapdos(counts_all = counts_all, counts_m = counts_m, n_resample = 100)
dt_unmeth <- zapdos(counts_all = counts_all, counts_m = counts_um, n_resample = 100)

### EXPORT 
fwrite(dt_meth, mae.meth.outfile)
fwrite(dt_unmeth, mae.unmeth.outfile)


#####

# caterpie <- function(sub){
# 
#   mutation = rep(sub$mutation[1], 2)
#   model = c("A-CG", "A-CG_meth")
#   MAE =  c(mean(sub$AE_null), mean(sub$AE))
#   MAE_CI95 = c(1.96 * (sd(sub$AE_null) / sqrt(length(sub$AE_null))),
#                1.96 * (sd(sub$AE) / sqrt(length(sub$AE))))
#   output <- data.table(mutation, model, MAE, MAE_CI95)
# 
#   return(output)
# }
# metapod_k1CGrel <- function(sub){
#   data.table(mutation = sub$mutation,
#              model = sub$model,
#              MAE_k1CG_rel = ((sub$MAE - sub$MAE[sub$model == "A-CG"]) / sub$MAE[sub$model == "A-CG"])*100,
#              MAE_k1CG_rel_CI95 = (sub$MAE_CI95 / sub$MAE[sub$model == "A-CG"])*100
#   )
# }
# 
# 
# dt_meth_list <- list()
# for (i in 1:length(unique(dt_meth$mutation))){
#   sub <- dt_meth[mutation == unique(dt_meth$mutation)[i]]
#   dt_meth_list[[i]] <- caterpie(sub)
# }
# dt_meth <- do.call("rbind", dt_meth_list)
# dt_meth[, 3:4] <- lapply(dt_meth[, 3:4], as.numeric) # convert
# for (i in 1:length(unique(dt_meth$mutation))){
#   sub <- dt_meth[mutation == unique(dt_meth$mutation)[i]]
#   dt_meth_list[[i]] <- metapod_k1CGrel(sub)
# }
# dt_meth <- do.call("rbind", dt_meth_list)
# dt_meth$model[dt_meth$model == "A-CG_meth"] <- "CG+"
# 
# dt_unmeth_list <- list()
# for (i in 1:length(unique(dt_unmeth$mutation))){
#   sub <- dt_unmeth[mutation == unique(dt_unmeth$mutation)[i]]
#   dt_unmeth_list[[i]] <- caterpie(sub)
# }
# dt_unmeth <- do.call("rbind", dt_unmeth_list)
# dt_unmeth[, 3:4] <- lapply(dt_unmeth[, 3:4], as.numeric) # convert to numeric
# for (i in 1:length(unique(dt_unmeth$mutation))){
#   sub <- dt_unmeth[mutation == unique(dt_unmeth$mutation)[i]]
#   dt_unmeth_list[[i]] <- metapod_k1CGrel(sub)
# }
# dt_unmeth <- do.call("rbind", dt_unmeth_list)
# dt_unmeth$model[dt_unmeth$model == "A-CG_meth"] <- "CG-"
# 
# dt_mae <- rbind(dt_meth, dt_unmeth)
# dt_mae <- dt_mae[model != "A-CG"]
# dt_mae$model <- factor(dt_mae$model, levels = c("CG+", "CG-"))
# 
# fig_mae <- ggplot(dt_mae, aes(x=model, y=MAE_k1CG_rel, group=mutation, color=mutation)) +
#   geom_errorbar(aes(ymax = MAE_k1CG_rel+MAE_k1CG_rel_CI95, ymin = MAE_k1CG_rel-MAE_k1CG_rel_CI95),
#                 stat = "identity",
#                 position = position_dodge(width = 0.5),
#                 width=0.2,
#                 size=1.2) +
#   geom_point(size = 1,
#              position = position_dodge(width = 0.5)) +
#   xlab("") +
#   ylab("MAE per kb\n(% change relative to null model)") +
#   # ggtitle() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 26, face = "bold"),
#         legend.position = c(0.8, 0.8),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 14),
#         legend.background = element_rect(linetype="solid", colour ="black"),
#         plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   )
# fig_mae







