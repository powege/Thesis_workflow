# for each substitution type (A>C, A>G, A>T, and C>A, C>G, C>T)
# plot distribution of mu rate
# plot kmer composition ~ mutation rate percentile
# test with poisson glm

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

### IMPORT 
ann_counts_chr <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/m.mus_GRC38_k7_counts_by_annotation_and_chr.csv.gz")
k7_psnv <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz")

### FORMAT

k7_psnv_any <- setDT(k7_psnv)[, .(k7_to_N = sum(k7_to_N)), by = .(k7_from, k7_from_N)]
k7_psnv_any$k7_mu_rate <- k7_psnv_any$k7_to_N / k7_psnv_any$k7_from_N
k7_psnv_any$to <- "X"
k7_psnv_any <- k7_psnv_any[,c("k7_from", "k7_mu_rate", "to")]
k7_psnv_specific <- k7_psnv[,c("k7_from", "to", "k7_mu_rate")]
k7_psnv <- rbind(k7_psnv_specific, k7_psnv_any)

ann_counts <- setDT(ann_counts_chr)[, .(category_k7_N = sum(k7_from_N)), by = .(k7_from, category)]
ann_counts <- ann_counts[category == "All"]
ann_counts$k7_percentage <- (ann_counts$category_k7_N / sum(ann_counts$category_k7_N))*100
dt <- k7_psnv[ann_counts[,c("k7_from", "k7_percentage", "category", "category_k7_N")], on = "k7_from"]
dt_nCG <- dt[which(stri_sub(dt$k7_from, 4, -3) != "CG")]
dt_CG <- dt[which(stri_sub(dt$k7_from, 4, -3) == "CG")] 
dt_nCG$mu_type <- paste0(stri_sub(dt_nCG$k7_from, 4, -4), ">", dt_nCG$to)
dt_CG$mu_type <- paste0(stri_sub(dt_CG$k7_from, 4, -4), "G>", dt_CG$to, "G")
dt <- rbind(dt_nCG, dt_CG)

# sub <- dt[mu_type == "A>C"]
frog <- function(sub){
  percentile <- ecdf(sub$k7_mu_rate)
  sub$k7_mu_rate_percentile <- ceiling(percentile(sub$k7_mu_rate)*100)
  return(sub)
}
dt <- ddply(dt, "mu_type", frog)
dt_ag <- aggregate(list(k7_N=dt$k7_percentage), by=list(percentile=dt$k7_mu_rate_percentile, mu_type=dt$mu_type), FUN=sum)


## HISTOGRAMS

hist(dt$k7_mu_rate[dt$mu_type == "A>X"])
hist(dt$k7_mu_rate[dt$mu_type == "A>C"])
hist(dt$k7_mu_rate[dt$mu_type == "A>G"])
hist(dt$k7_mu_rate[dt$mu_type == "A>T"])

hist(dt$k7_mu_rate[dt$mu_type == "C>X"])
hist(dt$k7_mu_rate[dt$mu_type == "C>A"])
hist(dt$k7_mu_rate[dt$mu_type == "C>G"])
hist(dt$k7_mu_rate[dt$mu_type == "C>T"])

hist(dt$k7_mu_rate[dt$mu_type == "CG>XG"])
hist(dt$k7_mu_rate[dt$mu_type == "CG>AG"])
hist(dt$k7_mu_rate[dt$mu_type == "CG>GG"])
hist(dt$k7_mu_rate[dt$mu_type == "CG>TG"])

## SCATTER

cor.test(dt_ag$k7_N[dt_ag$mu_type == "A>X"], dt_ag$percentile[dt_ag$mu_type == "A>X"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "A>C"], dt_ag$percentile[dt_ag$mu_type == "A>C"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "A>G"], dt_ag$percentile[dt_ag$mu_type == "A>G"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "A>T"], dt_ag$percentile[dt_ag$mu_type == "A>T"])

cor.test(dt_ag$k7_N[dt_ag$mu_type == "C>X"], dt_ag$percentile[dt_ag$mu_type == "C>X"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "C>A"], dt_ag$percentile[dt_ag$mu_type == "C>A"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "C>G"], dt_ag$percentile[dt_ag$mu_type == "C>G"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "C>T"], dt_ag$percentile[dt_ag$mu_type == "C>T"])

cor.test(dt_ag$k7_N[dt_ag$mu_type == "CG>XG"], dt_ag$percentile[dt_ag$mu_type == "CG>XG"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "CG>AG"], dt_ag$percentile[dt_ag$mu_type == "CG>AG"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "CG>GG"], dt_ag$percentile[dt_ag$mu_type == "CG>GG"])
cor.test(dt_ag$k7_N[dt_ag$mu_type == "CG>TG"], dt_ag$percentile[dt_ag$mu_type == "CG>TG"])

plot(dt_ag$k7_N[dt_ag$mu_type == "A>X"] ~ dt_ag$percentile[dt_ag$mu_type == "A>X"])
plot(dt_ag$k7_N[dt_ag$mu_type == "A>C"] ~ dt_ag$percentile[dt_ag$mu_type == "A>C"])
plot(dt_ag$k7_N[dt_ag$mu_type == "A>G"] ~ dt_ag$percentile[dt_ag$mu_type == "A>G"])
plot(dt_ag$k7_N[dt_ag$mu_type == "A>T"] ~ dt_ag$percentile[dt_ag$mu_type == "A>T"])

plot(dt_ag$k7_N[dt_ag$mu_type == "C>X"] ~ dt_ag$percentile[dt_ag$mu_type == "C>X"])
plot(dt_ag$k7_N[dt_ag$mu_type == "C>A"] ~ dt_ag$percentile[dt_ag$mu_type == "C>A"])
plot(dt_ag$k7_N[dt_ag$mu_type == "C>G"] ~ dt_ag$percentile[dt_ag$mu_type == "C>G"])
plot(dt_ag$k7_N[dt_ag$mu_type == "C>T"] ~ dt_ag$percentile[dt_ag$mu_type == "C>T"])

plot(dt_ag$k7_N[dt_ag$mu_type == "CG>XG"] ~ dt_ag$percentile[dt_ag$mu_type == "CG>XG"])
plot(dt_ag$k7_N[dt_ag$mu_type == "CG>AG"] ~ dt_ag$percentile[dt_ag$mu_type == "CG>AG"])
plot(dt_ag$k7_N[dt_ag$mu_type == "CG>GG"] ~ dt_ag$percentile[dt_ag$mu_type == "CG>GG"])
plot(dt_ag$k7_N[dt_ag$mu_type == "CG>TG"] ~ dt_ag$percentile[dt_ag$mu_type == "CG>TG"])

## POISSON GLM

summary(glm(dt$category_k7_N[dt$mu_type == "A>X"] ~ dt$k7_mu_rate[dt$mu_type == "A>X"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "A>C"] ~ dt$k7_mu_rate[dt$mu_type == "A>C"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "A>G"] ~ dt$k7_mu_rate[dt$mu_type == "A>G"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "A>T"] ~ dt$k7_mu_rate[dt$mu_type == "A>T"], family = "poisson"))

summary(glm(dt$category_k7_N[dt$mu_type == "C>X"] ~ dt$k7_mu_rate[dt$mu_type == "C>X"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "C>A"] ~ dt$k7_mu_rate[dt$mu_type == "C>A"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "C>G"] ~ dt$k7_mu_rate[dt$mu_type == "C>G"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "C>T"] ~ dt$k7_mu_rate[dt$mu_type == "C>T"], family = "poisson"))

summary(glm(dt$category_k7_N[dt$mu_type == "CG>XG"] ~ dt$k7_mu_rate[dt$mu_type == "CG>XG"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "CG>AG"] ~ dt$k7_mu_rate[dt$mu_type == "CG>AG"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "CG>GG"] ~ dt$k7_mu_rate[dt$mu_type == "CG>GG"], family = "poisson"))
summary(glm(dt$category_k7_N[dt$mu_type == "CG>TG"] ~ dt$k7_mu_rate[dt$mu_type == "CG>TG"], family = "poisson"))



######

# rm(list = ls())
# graphics.off()
# 
# library(data.table)
# 
# ### IMPORT 
# ann_counts_chr <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/m.mus_GRC38_k7_counts_by_annotation_and_chr.csv.gz")
# k7_psnv <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz")
# 
# ### FORMAT
# 
# k7_psnv <- setDT(k7_psnv)[, .(k7_to_N = sum(k7_to_N)), by = .(k7_from, k7_from_N)]
# k7_psnv$k7_mu_rate <- k7_psnv$k7_to_N / k7_psnv$k7_from_N
# 
# ann_counts <- setDT(ann_counts_chr)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]
# ann_counts_wide <- dcast(ann_counts, k7_from ~ category, value.var = "k7_from_N")
# percentage_function <- function(vec){( vec/sum(vec) )*100}
# ann_counts_wide2 <- cbind(ann_counts_wide[,1], apply(ann_counts_wide[,2:ncol(ann_counts_wide)], 2, percentage_function))
# ann_counts2 <- melt(ann_counts_wide2, id.vars = c("k7_from"))
# 
# ### DISTRIBUTION
# 
# # hist(k7_psnv$k7_mu_rate)
# A_psnv <- k7_psnv[which(stri_sub(k7_psnv$k7_from, 4, -4) == "A")]
# hist(A_psnv$k7_mu_rate)
# C_psnv <- k7_psnv[which(stri_sub(k7_psnv$k7_from, 4, -4) == "C")]
# hist(C_psnv$k7_mu_rate)
# C_nonCG_psnv <- C_psnv[which(stri_sub(C_psnv$k7_from, 4, -3) != "CG")]
# hist(C_nonCG_psnv$k7_mu_rate)
# C_CG_psnv <- k7_psnv[which(stri_sub(k7_psnv$k7_from, 4, -3) == "CG")]
# hist(C_CG_psnv$k7_mu_rate)
# 
# ann_counts_all_A <- ann_counts[category == "All"]
# ann_counts_all_A <- ann_counts_all_A[which(stri_sub(ann_counts_all_A$k7_from, 4, -4) == "A")]
# sum(ann_counts_all_A$k7_from_N)
# ann_counts_all_C <- ann_counts[category == "All"]
# ann_counts_all_C <- ann_counts_all_C[which(stri_sub(ann_counts_all_C$k7_from, 4, -4) == "C")]
# sum(ann_counts_all_C$k7_from_N)
# 
# ### CORRELATION
# 
# # test correlations
# dt_AC_wide <- ann_counts_wide[k7_psnv[,c("k7_from", "k7_mu_rate")], on = c("k7_from")]
# cor(dt_AC_wide[,2:12])
# dt_nonCG_wide <- dt_AC_wide[which(stri_sub(dt_AC_wide$k7_from, 4, -3) != "CG")]
# cor(dt_nonCG_wide[,2:12])
# dt_A_wide <- dt_AC_wide[which(stri_sub(dt_AC_wide$k7_from, 4, -4) == "A"),]
# cor(dt_A_wide[,2:12])
# dt_C_nonCG_wide <- dt_AC_wide[which(stri_sub(dt_AC_wide$k7_from, 4, -4) == "C"),]
# dt_C_nonCG_wide <- dt_C_nonCG_wide[which(stri_sub(dt_C_nonCG_wide$k7_from, 4, -3) != "CG"),]
# cor(dt_C_nonCG_wide[,2:12])
# dt_CG_wide <- dt_AC_wide[which(stri_sub(dt_AC_wide$k7_from, 4, -3) == "CG")]
# cor(dt_CG_wide[,2:12])
# 
# # plot(dt_A_wide$All ~ dt_A_wide$k7_mu_rate_percentile)
# # plot(dt_C_nonCG_wide$All ~ dt_C_nonCG_wide$k7_mu_rate_percentile)
# # plot(dt_CG_wide$All ~ dt_CG_wide$k7_mu_rate_percentile)
# # plot(dt_A_wide$All ~ dt_A_wide$k7_mu_rate)
# # plot(dt_C_nonCG_wide$All ~ dt_C_nonCG_wide$k7_mu_rate)
# # plot(dt_CG_wide$All ~ dt_CG_wide$k7_mu_rate_percentile)
# 
# # glm
# # summary(glm(formula = All ~ k7_mu_rate, data = dt_nonCG_wide, family = "poisson"))
# summary(glm(formula = All ~ k7_mu_rate, data = dt_CG_wide, family = "poisson"))
# summary(glm(formula = All ~ k7_mu_rate, data = dt_A_wide, family = "poisson"))
# summary(glm(formula = All ~ k7_mu_rate, data = dt_C_nonCG_wide, family = "poisson"))
# 
# dt_AC <- ann_counts2[k7_psnv[,c("k7_from", "k7_mu_rate")], on = c("k7_from")]
# dt_nonCG <- dt_AC[which(stri_sub(dt_AC$k7_from, 4, -3) != "CG")]
# dt_A <- dt_AC[which(stri_sub(dt_AC$k7_from, 4, -4) == "A"),]
# dt_C_nonCG <- dt_AC[which(stri_sub(dt_AC$k7_from, 4, -4) == "C"),]
# dt_C_nonCG <- dt_C_nonCG[which(stri_sub(dt_C_nonCG$k7_from, 4, -3) != "CG"),]
# dt_CG <- dt_AC[which(stri_sub(dt_AC$k7_from, 4, -3) == "CG")]
# 
# # ggplot(dt_A, aes(x=k7_mu_rate, y=value, color=variable)) +
# #   geom_point(alpha = 1/10) + 
# #   geom_smooth(method=lm)
# # 
# # ggplot(dt_C_nonCG, aes(x=k7_mu_rate, y=value, color=variable)) +
# #   geom_point(alpha = 1/10) + 
# #   geom_smooth(method=lm)
# # 
# # ggplot(dt_CG, aes(x=k7_mu_rate, y=value, color=variable)) +
# #   geom_point(alpha = 1/10) + 
# #   geom_smooth(method=lm)
# 
# ### OR 
# 
# # dt_in <- dt_A_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")]
# # percentiles <- 1:10
# # n_in_group <- nrow(dt_A_wide[dt_A_wide$k7_mu_rate_percentile %in% percentiles])
# # n_out_group <- nrow(dt_A_wide[!dt_A_wide$k7_mu_rate_percentile %in% percentiles])
# 
# odds_ratio <- function(dt_in, percentiles, n_in_group, n_out_group){
#   
#   groups <- unique(dt_in$category)
#   colnames(dt_in) <- c("k7_from_N", "category", "percentile")
#   percentiles_vec <- rep(paste(min(percentiles), "to", max(percentiles), sep = " "), length(groups))
#   A <- rep(NA, length(groups))
#   B <- rep(NA, length(groups))
#   C <- rep(NA, length(groups))
#   D <- rep(NA, length(groups))
#   
#   for (i in 1:length(groups)){
#     A[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i] &
#                                   dt_in$percentile %in% percentiles])
#     B[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i]]) * (n_in_group / (n_in_group + n_out_group))
#     C[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i] &
#                                   !dt_in$percentile %in% percentiles])
#     D[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i]]) * (n_out_group / (n_in_group + n_out_group))
#   }
#   OR <- (A/B)/(C/D)
#   CI95 <- 1.96*sqrt((1/A) + (1/(B)) + (1/C) + (1/D))
#   out <- data.table(percentiles = percentiles_vec, 
#                     category = groups, 
#                     A, B, C, D, OR, CI95)
#   
#   return(out)
# }
# 
# percentile <- ecdf(dt_nonCG_wide$k7_mu_rate)
# dt_nonCG_wide$k7_mu_rate_percentile <- ceiling(percentile(dt_nonCG_wide$k7_mu_rate)*100)
# dt_nonCG_long <- melt(dt_nonCG_wide,
#                       id.vars=c("k7_mu_rate", "k7_mu_rate_percentile"),
#                       measure.vars=c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter"),
#                       variable.name="category",
#                       value.name="k7_from_N"
# )
# or_nonCG_low <- odds_ratio(dt_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 1:10, nrow(dt_nonCG_wide[dt_nonCG_wide$k7_mu_rate_percentile %in% 1:10]), nrow(dt_nonCG_wide[!dt_nonCG_wide$k7_mu_rate_percentile %in% 1:10]))
# or_nonCG_high <- odds_ratio(dt_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 91:100, nrow(dt_nonCG_wide[dt_nonCG_wide$k7_mu_rate_percentile %in% 91:100]), nrow(dt_nonCG_wide[!dt_nonCG_wide$k7_mu_rate_percentile %in% 91:100]))
# or_nonCG_mid <- odds_ratio(dt_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 11:90, nrow(dt_nonCG_wide[dt_nonCG_wide$k7_mu_rate_percentile %in% 11:90]), nrow(dt_nonCG_wide[!dt_nonCG_wide$k7_mu_rate_percentile %in% 11:90]))
# 
# percentile <- ecdf(dt_CG_wide$k7_mu_rate)
# dt_CG_wide$k7_mu_rate_percentile <- ceiling(percentile(dt_CG_wide$k7_mu_rate)*100)
# dt_CG_long <- melt(dt_CG_wide,
#                    id.vars=c("k7_mu_rate", "k7_mu_rate_percentile"),
#                    measure.vars=c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter"),
#                    variable.name="category",
#                    value.name="k7_from_N"
# )
# or_CG_low <- odds_ratio(dt_CG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 1:10, nrow(dt_CG_wide[dt_CG_wide$k7_mu_rate_percentile %in% 1:10]), nrow(dt_CG_wide[!dt_CG_wide$k7_mu_rate_percentile %in% 1:10]))
# or_CG_high <- odds_ratio(dt_CG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 91:100, nrow(dt_CG_wide[dt_CG_wide$k7_mu_rate_percentile %in% 91:100]), nrow(dt_CG_wide[!dt_CG_wide$k7_mu_rate_percentile %in% 91:100]))
# or_CG_mid <- odds_ratio(dt_CG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 11:90, nrow(dt_CG_wide[dt_CG_wide$k7_mu_rate_percentile %in% 11:90]), nrow(dt_CG_wide[!dt_CG_wide$k7_mu_rate_percentile %in% 11:90]))
# 
# percentile <- ecdf(dt_A_wide$k7_mu_rate)
# dt_A_wide$k7_mu_rate_percentile <- ceiling(percentile(dt_A_wide$k7_mu_rate)*100)
# dt_A_long <- melt(dt_A_wide, 
#                   id.vars=c("k7_mu_rate", "k7_mu_rate_percentile"),
#                   measure.vars=c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter"),
#                   variable.name="category",
#                   value.name="k7_from_N"
# )
# or_A_low <- odds_ratio(dt_A_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 1:10, nrow(dt_A_wide[dt_A_wide$k7_mu_rate_percentile %in% 1:10]), nrow(dt_A_wide[!dt_A_wide$k7_mu_rate_percentile %in% 1:10]))
# or_A_high <- odds_ratio(dt_A_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 91:100, nrow(dt_A_wide[dt_A_wide$k7_mu_rate_percentile %in% 91:100]), nrow(dt_A_wide[!dt_A_wide$k7_mu_rate_percentile %in% 91:100]))
# or_A_mid <- odds_ratio(dt_A_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 11:90, nrow(dt_A_wide[dt_A_wide$k7_mu_rate_percentile %in% 11:90]), nrow(dt_A_wide[!dt_A_wide$k7_mu_rate_percentile %in% 11:90]))
# 
# percentile <- ecdf(dt_C_nonCG_wide$k7_mu_rate)
# dt_C_nonCG_wide$k7_mu_rate_percentile <- ceiling(percentile(dt_C_nonCG_wide$k7_mu_rate)*100)
# dt_C_nonCG_long <- melt(dt_C_nonCG_wide, 
#                         id.vars=c("k7_mu_rate", "k7_mu_rate_percentile"),
#                         measure.vars=c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter"),
#                         variable.name="category",
#                         value.name="k7_from_N"
# )
# or_C_nonCG_low <- odds_ratio(dt_C_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 1:10, nrow(dt_C_nonCG_wide[dt_C_nonCG_wide$k7_mu_rate_percentile %in% 1:10]), nrow(dt_C_nonCG_wide[!dt_C_nonCG_wide$k7_mu_rate_percentile %in% 1:10]))
# or_C_nonCG_high <- odds_ratio(dt_C_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 91:100, nrow(dt_C_nonCG_wide[dt_C_nonCG_wide$k7_mu_rate_percentile %in% 91:100]), nrow(dt_C_nonCG_wide[!dt_C_nonCG_wide$k7_mu_rate_percentile %in% 91:100]))
# or_C_nonCG_mid <- odds_ratio(dt_C_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 11:90, nrow(dt_C_nonCG_wide[dt_C_nonCG_wide$k7_mu_rate_percentile %in% 11:90]), nrow(dt_C_nonCG_wide[!dt_C_nonCG_wide$k7_mu_rate_percentile %in% 11:90]))
# 
# 
# ### BOXPLOT 
# 
# # boxplot mutation rates (all)
# # A > 
# # C (nonCG) >
# # C (CG) >
# # All >
# 
# # boxplot mutation rates (specific)
# # A > C
# # C (nonCG) > A
# # C (CG) > A
# # A > G
# # C (nonCG) > G
# # C (CG) > G
# # A > T
# # C (nonCG) > T
# # C (CG) > T
# 
# #####
# categories <- c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter")
# i = 1
# dt_AC_all <- dt_AC[variable == categories[i]]
# percentile <- ecdf(dt_AC_all$k7_mu_rate)
# dt_AC_all$k7_mu_rate_percentile <- ceiling(percentile(dt_AC_all$k7_mu_rate)*100)
# dt_AC_all$from <- stri_sub(dt_AC_all$k7_from, 4, -4)
# dt_AC_all_ag <- aggregate(list(kmer_percentage=dt_AC_all$value), by=list(percentile=dt_AC_all$k7_mu_rate_percentile), FUN=sum)
# plot(dt_AC_all_ag$kmer_percentage ~ dt_AC_all_ag$percentile)
# 
# tmp_1 <- dt_AC_all[k7_mu_rate_percentile <= 48]
# table(tmp_1$from)
# tmp_2 <- dt_AC_all[k7_mu_rate_percentile >= 49 & k7_mu_rate_percentile <= 87]
# table(tmp_2$from)
# 
# dt_nonCG_all <- dt_nonCG_long[category == categories[i]]
# dt_nonCG_all <- aggregate(list(n=dt_nonCG_all$k7_from_N), by=list(percentile=dt_nonCG_all$k7_mu_rate_percentile), FUN=sum)
# dt_nonCG_all$kmer_percentage <- (dt_nonCG_all$n / sum(dt_nonCG_all$n)) * 100
# plot(dt_nonCG_all$kmer_percentage ~ dt_nonCG_all$percentile)
# 
# dt_A_all <- dt_A_long[category == categories[i]]
# dt_A_all <- aggregate(list(n=dt_A_all$k7_from_N), by=list(percentile=dt_A_all$k7_mu_rate_percentile), FUN=sum)
# dt_A_all$kmer_percentage <- (dt_A_all$n / sum(dt_A_all$n)) * 100
# plot(dt_A_all$kmer_percentage ~ dt_A_all$percentile)
# 
# dt_C_nonCG_all <- dt_C_nonCG_long[category == categories[i]]
# dt_C_nonCG_all <- aggregate(list(n=dt_C_nonCG_all$k7_from_N), by=list(percentile=dt_C_nonCG_all$k7_mu_rate_percentile), FUN=sum)
# dt_C_nonCG_all$kmer_percentage <- (dt_C_nonCG_all$n / sum(dt_C_nonCG_all$n)) * 100
# plot(dt_C_nonCG_all$kmer_percentage ~ dt_C_nonCG_all$percentile)
# 
# dt_C_CG_all <- dt_CG_long[category == categories[i]]
# dt_C_CG_all <- aggregate(list(n=dt_C_CG_all$k7_from_N), by=list(percentile=dt_C_CG_all$k7_mu_rate_percentile), FUN=sum)
# dt_C_CG_all$kmer_percentage <- (dt_C_CG_all$n / sum(dt_C_CG_all$n)) * 100
# plot(dt_C_CG_all$kmer_percentage ~ dt_C_CG_all$percentile)


#####




