rm(list = ls())
graphics.off()

library(data.table)

### IMPORT
m_ann_counts_chr <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/m.mus_GRC38_k7_counts_by_annotation_and_chr.csv.gz")
m_k7_psnv <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz")
h_ann_counts_chr <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/h.sap_GRC38_k7_counts_by_annotation_and_chr.csv.gz")
h_k7_psnv <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/Mouse_SNV_rates/Aggarwala_Voight_2016_7mer_pmu.csv")

### FORMAT

# mouse pSNV
m_k7_psnv <- setDT(m_k7_psnv)[, .(k7_to_N = sum(k7_to_N)), by = .(k7_from, k7_from_N)]
m_k7_psnv$mmus_mu_rate <- m_k7_psnv$k7_to_N / m_k7_psnv$k7_from_N
m_k7_psnv <- m_k7_psnv[,c("k7_from", "mmus_mu_rate")]

# human pSNV
colnames(h_k7_psnv) <- c("k7_from", "k7_to", "African_pmu", "Asian_pmu", "European_pmu", "k7_from_reverse", "k7_to_reverse")
h_k7_psnv <- h_k7_psnv[,c("k7_from", "African_pmu", "Asian_pmu", "European_pmu")]
h_k7_psnv$k7_from_N <- 100
h_k7_psnv$k7_to_N <- h_k7_psnv$African_pmu * h_k7_psnv$k7_from_N
h_k7_psnv <- setDT(h_k7_psnv)[, .(k7_to_N = sum(k7_to_N)), by = .(k7_from, k7_from_N)]
h_k7_psnv$hsap_mu_rate <- h_k7_psnv$k7_to_N / h_k7_psnv$k7_from_N
h_k7_psnv <- h_k7_psnv[,c("k7_from", "hsap_mu_rate")]

k7_psnv <- h_k7_psnv[m_k7_psnv, on = "k7_from"]

# k7 counts by ann
h_ann_counts <- setDT(h_ann_counts_chr)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]
m_ann_counts <- setDT(m_ann_counts_chr)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]
colnames(m_ann_counts) <- c("k7_from", "category", "mmus_from_N")
colnames(h_ann_counts) <- c("k7_from", "category", "hsap_from_N")
ann_counts <- m_ann_counts[h_ann_counts, on = c("k7_from", "category")]
ann_counts <- ann_counts[category == "All"]
ann_counts$mmus_k7_percentage <- (ann_counts$mmus_from_N / sum(ann_counts$mmus_from_N))*100
ann_counts$hsap_k7_percentage <- (ann_counts$hsap_from_N / sum(ann_counts$hsap_from_N))*100
# plot(ann_counts$mmus_k7_percentage, ann_counts$hsap_k7_percentage)

dt <- ann_counts[k7_psnv, on = "k7_from"]
dt_nCG <- dt[which(stri_sub(dt$k7_from, 4, -3) != "CG")]
dt_CG <- dt[which(stri_sub(dt$k7_from, 4, -3) == "CG")] 
dt_nCG$mu_type <- paste0(stri_sub(dt_nCG$k7_from, 4, -4), ">", dt_nCG$to)
dt_CG$mu_type <- paste0(stri_sub(dt_CG$k7_from, 4, -4), "G>", dt_CG$to, "")
dt <- rbind(dt_nCG, dt_CG)

percentile <- ecdf(dt$hsap_mu_rate)
dt$hsap_mu_rate_percentile <- percentile(dt$hsap_mu_rate)
percentile <- ecdf(dt$mmus_mu_rate)
dt$mmus_mu_rate_percentile <- percentile(dt$mmus_mu_rate)
percentile <- ecdf(dt$hsap_k7_percentage)
dt$hsap_k7_percentage_percentile <- percentile(dt$hsap_k7_percentage)
percentile <- ecdf(dt$mmus_k7_percentage)
dt$mmus_k7_percentage_percentile <- percentile(dt$mmus_k7_percentage)

dt$k7_percentage_dif <- dt$mmus_k7_percentage - dt$hsap_k7_percentage
dt$mu_rate_dif <- dt$mmus_mu_rate - dt$hsap_mu_rate
dt$k7_percentage_pdif <- dt$mmus_k7_percentage_percentile - dt$hsap_k7_percentage_percentile
dt$mu_rate_pdif <- dt$mmus_mu_rate_percentile - dt$hsap_mu_rate_percentile

plot(dt$mu_rate_dif[dt$mu_type == "A>"], dt$k7_percentage_dif[dt$mu_type == "A>"])
plot(dt$mu_rate_dif[dt$mu_type == "C>"], dt$k7_percentage_dif[dt$mu_type == "C>"])
plot(dt$mu_rate_dif[dt$mu_type == "CG>"], dt$k7_percentage_dif[dt$mu_type == "CG>"])

cor.test(dt$mu_rate_dif[dt$mu_type == "A>"], dt$k7_percentage_dif[dt$mu_type == "A>"])
cor.test(dt$mu_rate_dif[dt$mu_type == "C>"], dt$k7_percentage_dif[dt$mu_type == "C>"])
cor.test(dt$mu_rate_dif[dt$mu_type == "CG>"], dt$k7_percentage_dif[dt$mu_type == "CG>"])

plot(dt$mu_rate_pdif[dt$mu_type == "A>"], dt$k7_percentage_pdif[dt$mu_type == "A>"])
plot(dt$mu_rate_pdif[dt$mu_type == "C>"], dt$k7_percentage_pdif[dt$mu_type == "C>"])
plot(dt$mu_rate_pdif[dt$mu_type == "CG>"], dt$k7_percentage_pdif[dt$mu_type == "CG>"])

cor.test(dt$mu_rate_pdif[dt$mu_type == "A>"], dt$k7_percentage_pdif[dt$mu_type == "A>"])
cor.test(dt$mu_rate_pdif[dt$mu_type == "C>"], dt$k7_percentage_pdif[dt$mu_type == "C>"])
cor.test(dt$mu_rate_pdif[dt$mu_type == "CG>"], dt$k7_percentage_pdif[dt$mu_type == "CG>"])


#####


# k7_psnv_nonCG <- k7_psnv[which(stri_sub(k7_psnv$k7_from, 4, -3) != "CG"),]
# k7_psnv_nonCG <- k7_psnv[which(stri_sub(k7_psnv$k7_from, 4, -4) == "C" & stri_sub(k7_psnv$k7_from, 4, -3) != "CG"),]
k7_psnv_nonCG <- k7_psnv[which(stri_sub(k7_psnv$k7_from, 4, -4) == "A"),]
k7_psnv_nonCG$hsap_mu_rate_rel <- (k7_psnv_nonCG$hsap_mu_rate - mean(k7_psnv_nonCG$hsap_mu_rate)) / mean(k7_psnv_nonCG$hsap_mu_rate)
k7_psnv_nonCG$mmus_mu_rate_rel <- (k7_psnv_nonCG$mmus_mu_rate - mean(k7_psnv_nonCG$mmus_mu_rate)) / mean(k7_psnv_nonCG$mmus_mu_rate)
# plot(k7_psnv_nonCG$hsap_mu_rate_rel, k7_psnv_nonCG$mmus_mu_rate_rel)
k7_psnv_nonCG$mmus_mu_rate_dif <- k7_psnv_nonCG$hsap_mu_rate_rel - k7_psnv_nonCG$mmus_mu_rate_rel
# hist(k7_psnv_nonCG$mmus_mu_rate_dif)
# hist(k7_psnv_nonCG$mmus_mu_rate_dif[k7_psnv_nonCG$mmus_mu_rate_dif > -1000 & k7_psnv_nonCG$mmus_mu_rate_dif < 1000])
percentile <- ecdf(k7_psnv_nonCG$mmus_mu_rate_dif)
k7_psnv_nonCG$percentile <- ceiling(percentile(k7_psnv_nonCG$mmus_mu_rate_dif)*100)

# k7 counts by ann
h_ann_counts <- setDT(h_ann_counts_chr)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]
m_ann_counts <- setDT(m_ann_counts_chr)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]
colnames(m_ann_counts) <- c("k7_from", "category", "mmus_from_N")
colnames(h_ann_counts) <- c("k7_from", "category", "hsap_from_N")
ann_counts <- m_ann_counts[h_ann_counts, on = c("k7_from", "category")]

# ann_counts_nonCG <- ann_counts[which(stri_sub(ann_counts$k7_from, 4, -3) != "CG"),]
# ann_counts_nonCG <- ann_counts[which(stri_sub(ann_counts$k7_from, 4, -4) == "C" & stri_sub(ann_counts$k7_from, 4, -3) != "CG"),]
ann_counts_nonCG <- ann_counts[which(stri_sub(ann_counts$k7_from, 4, -4) == "A"),]
ann_counts_nonCG$mmus_from_N_rel <- (ann_counts_nonCG$mmus_from_N / sum(ann_counts_nonCG$mmus_from_N))*100
ann_counts_nonCG$hsap_from_N_rel <- (ann_counts_nonCG$hsap_from_N / sum(ann_counts_nonCG$hsap_from_N))*100
# plot(ann_counts_nonCG$hsap_from_N_rel, ann_counts_nonCG$mmus_from_N_rel)
ann_counts_nonCG$k7_from_N_dif <- ann_counts_nonCG$hsap_from_N_rel - ann_counts_nonCG$mmus_from_N_rel
# hist(ann_counts_nonCG$k7_from_N_dif)

tmp <- ann_counts_nonCG[category == "Enhancer - distal", c("k7_from", "k7_from_N_dif")][k7_psnv_nonCG[,c("k7_from", "mmus_mu_rate_dif")], on = "k7_from"]
cor.test(tmp$k7_from_N_dif, tmp$mmus_mu_rate_dif)
plot(tmp$k7_from_N_dif, tmp$mmus_mu_rate_dif)


for (i in 1:length(unique(ann_counts_nonCG$category))){
  ann_sub <- ann_counts_nonCG[category == unique(ann_counts_nonCG$category)[i]]

  k7_top <- k7_psnv_nonCG$k7_from[k7_psnv_nonCG$percentile >= 91] # greater in mouse
  k7_bottom <- k7_psnv_nonCG$k7_from[k7_psnv_nonCG$percentile <= 10] # greater in human
  
  A <- sum(ann_sub$mmus_from_N[ann_sub$k7_from %in% k7_bottom])
  B <- sum(ann_sub$mmus_from_N)
  C <- sum(ann_sub$hsap_from_N[ann_sub$k7_from %in% k7_bottom])
  D <- sum(ann_sub$hsap_from_N)
  
  # A <- sum(ann_sub$mmus_from_N[ann_sub$k7_from %in% k7_top])
  # B <- sum(ann_sub$mmus_from_N)
  # C <- sum(ann_sub$hsap_from_N[ann_sub$k7_from %in% k7_top])
  # D <- sum(ann_sub$hsap_from_N)
  
  OR <- (A/B)/(C/D)
  CI95 <- 1.96*sqrt((1/A) + (1/(B)) + (1/C) + (1/D))
}




#####

h_ann_counts <- h_ann_counts[category == "All"][,c("k7_from", "k7_from_N")]
m_ann_counts <- m_ann_counts[category == "All"][,c("k7_from", "k7_from_N")]
colnames(m_ann_counts) <- c("k7_from", "mmus_from_N")
colnames(h_ann_counts) <- c("k7_from", "hsap_from_N")
ann_counts <- m_ann_counts[h_ann_counts, on = "k7_from"]
ann_counts_nonCG <- ann_counts[which(stri_sub(ann_counts$k7_from, 4, -3) != "CG"),]
ann_counts_nonCG$mmus_from_N_rel <- (ann_counts_nonCG$mmus_from_N / sum(ann_counts_nonCG$mmus_from_N))*100
ann_counts_nonCG$hsap_from_N_rel <- (ann_counts_nonCG$hsap_from_N / sum(ann_counts_nonCG$hsap_from_N))*100
# plot(ann_counts_nonCG$hsap_from_N_rel, ann_counts_nonCG$mmus_from_N_rel)
ann_counts_nonCG$mmus_from_N_pdif <- ( (ann_counts_nonCG$hsap_from_N_rel - ann_counts_nonCG$mmus_from_N_rel) /  ann_counts_nonCG$mmus_from_N_rel ) *100
hist(ann_counts_nonCG$mmus_from_N_pdif)

tmp <- ann_counts_nonCG[,c("k7_from", "mmus_from_N_pdif")][k7_psnv_nonCG[,c("k7_from", "mmus_mu_rate_pdif")], on = "k7_from"]
cor.test(tmp$mmus_from_N_pdif, tmp$mmus_mu_rate_pdif, method = "spearman")
plot(tmp$mmus_from_N_pdif, tmp$mmus_mu_rate_pdif)
tmp2 <- tmp[mmus_mu_rate_pdif > -1000 & mmus_mu_rate_pdif < 1000,]
cor.test(tmp2$mmus_from_N_pdif, tmp2$mmus_mu_rate_pdif, method = "spearman")
plot(tmp2$mmus_from_N_pdif, tmp2$mmus_mu_rate_pdif)
#####


plot(k7_psnv_nonCG$hsap_mu_rate, k7_psnv_nonCG$mmus_mu_rate)
# plot(k7_psnv$hsap_mu_rate_rel, k7_psnv$mmus_mu_rate_rel)
cor.test(k7_psnv$hsap_mu_rate, k7_psnv$mmus_mu_rate)

plot(k7_psnv$hsap_mu_rate[which(stri_sub(k7_psnv$k7_from, 4, -3) != "CG")], k7_psnv$mmus_mu_rate[which(stri_sub(k7_psnv$k7_from, 4, -3) != "CG")])
cor.test(k7_psnv$hsap_mu_rate[which(stri_sub(k7_psnv$k7_from, 4, -3) != "CG")], k7_psnv$mmus_mu_rate[which(stri_sub(k7_psnv$k7_from, 4, -3) != "CG")])

h_ann_counts <- setDT(h_ann_counts_chr)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]
m_ann_counts <- setDT(m_ann_counts_chr)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]

h_ann_counts_wide <- dcast(h_ann_counts, k7_from ~ category, value.var = "k7_from_N")
m_ann_counts_wide <- dcast(m_ann_counts, k7_from ~ category, value.var = "k7_from_N")

percent_function <- function(vec){vec/sum(vec)*100}
h_ann_counts_wide_rel <- cbind(h_ann_counts_wide[,1], apply(h_ann_counts_wide[,2:ncol(h_ann_counts_wide)], 2, percent_function))
m_ann_counts_wide_rel <- cbind(m_ann_counts_wide[,1], apply(m_ann_counts_wide[,2:ncol(m_ann_counts_wide)], 2, percent_function))





#####

ann_counts <- setDT(ann_counts_chr)[, .(k7_from_N = sum(k7_from_N)), by = .(k7_from, category)]
ann_counts_wide <- dcast(ann_counts, k7_from ~ category, value.var = "k7_from_N")

dt_AC <- ann_counts_wide[k7_psnv[,c("k7_from", "k7_mu_rate")], on = c("k7_from")]
cor(dt_AC[,2:12])

dt_CG <- dt_AC[which(stri_sub(dt_AC$k7_from, 4, -3) == "CG")] 
cor(dt_CG[,2:12])
dt_nonCG <- dt_AC[which(stri_sub(dt_AC$k7_from, 4, -3) != "CG")] 
cor(dt_nonCG[,2:12])

dt_A <- dt_AC[which(stri_sub(dt_AC$k7_from, 4, -4) == "A"),]
cor(dt_A[,2:12])
dt_C_nonCG <- dt_nonCG[which(stri_sub(dt_nonCG$k7_from, 4, -4) == "C"),]
cor(dt_C_nonCG[,2:12])

plot(dt_A$k7_mu_rate, dt_A$All)
plot(dt_A$k7_mu_rate, dt_A$GERP)

plot(dt_C_nonCG$k7_mu_rate, dt_C_nonCG$All)
plot(dt_C_nonCG$k7_mu_rate, dt_C_nonCG$GERP)

plot(dt_CG$k7_mu_rate, dt_CG$All)
plot(dt_CG$k7_mu_rate, dt_CG$GERP)

plot(dt_nonCG$k7_mu_rate, dt_nonCG$All)
plot(dt_nonCG$k7_mu_rate, dt_nonCG$GERP)

hist(dt_A$k7_mu_rate)
hist(dt_A$All)
hist(dt_C_nonCG$k7_mu_rate)
hist(dt_C_nonCG$All)

summary(glm(formula = GERP ~ k7_mu_rate, data = dt_nonCG, family = "poisson"))
summary(glm(formula = GERP ~ k7_mu_rate, data = dt_CG, family = "poisson"))
summary(glm(formula = GERP ~ k7_mu_rate, data = dt_A, family = "poisson"))
summary(glm(formula = GERP ~ k7_mu_rate, data = dt_C_nonCG, family = "poisson"))


### OR 

# dt_in <- dt_A_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")]
# threshold <- 10
# highlow <- "<="
# in_group_n <- nrow(dt_A[dt_A$k7_mu_rate_percentile <= 10])
# out_group_n <- nrow(dt_A[dt_A$k7_mu_rate_percentile > 10])

odds_ratio <- function(dt_in, threshold, highlow, in_group_n, out_group_n){
  
  groups <- unique(dt_in$category)
  colnames(dt_in) <- c("k7_from_N", "category", "percentile")
  limit <- rep(threshold, length(groups))
  A <- rep(NA, length(groups))
  B <- rep(NA, length(groups))
  C <- rep(NA, length(groups))
  D <- rep(NA, length(groups))

  for (i in 1:length(groups)){
    if (highlow == ">="){
      A[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i] &
                                                       dt_in$percentile >= threshold])
      B[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i]]) * (in_group_n / (in_group_n + out_group_n))
      C[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i] &
                                                       dt_in$percentile < threshold])
      D[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i]]) * (out_group_n / (in_group_n + out_group_n))
    }
    if (highlow == "<="){
      
      A[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i] &
                                    dt_in$percentile <= threshold])
      B[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i]]) * (in_group_n / (in_group_n + out_group_n))
      C[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i] &
                                    dt_in$percentile > threshold])
      D[i] <- sum(dt_in$k7_from_N[dt_in$category %in% groups[i]]) * (out_group_n / (in_group_n + out_group_n))
    }
  }
  OR <- (A/B)/(C/D)
  CI95 <- 1.96*sqrt((1/A) + (1/(B)) + (1/C) + (1/D))
  out <- data.table(limit, groups, A, B, C, D, OR, CI95)
  
  return(out)
}

percentile <- ecdf(dt_nonCG$k7_mu_rate)
dt_nonCG$k7_mu_rate_percentile <- ceiling(percentile(dt_nonCG$k7_mu_rate)*100)
dt_nonCG_long <- melt(dt_nonCG,
                      id.vars=c("k7_mu_rate", "k7_mu_rate_percentile"),
                      measure.vars=c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter"),
                      variable.name="category",
                      value.name="k7_from_N"
)
or_nonCG_low <- odds_ratio(dt_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 10, "<=", nrow(dt_nonCG[dt_nonCG$k7_mu_rate_percentile <= 10]), nrow(dt_nonCG[dt_nonCG$k7_mu_rate_percentile > 10]))
or_nonCG_high <- odds_ratio(dt_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 91, ">=", nrow(dt_nonCG[dt_nonCG$k7_mu_rate_percentile >= 91]), nrow(dt_nonCG[dt_nonCG$k7_mu_rate_percentile < 91]))

# percentile <- ecdf(dt_CG$k7_mu_rate)
# dt_CG$k7_mu_rate_percentile <- ceiling(percentile(dt_CG$k7_mu_rate)*100)
# dt_CG_long <- melt(dt_CG, 
#                    id.vars=c("k7_mu_rate", "k7_mu_rate_percentile"),
#                    measure.vars=c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter"),
#                    variable.name="category",
#                    value.name="k7_from_N"
# )
# or_CG_low <- odds_ratio(dt_CG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 10, "<=", nrow(dt_CG[dt_CG$k7_mu_rate_percentile <= 10]), nrow(dt_CG[dt_CG$k7_mu_rate_percentile > 10]))
# or_CG_high <- odds_ratio(dt_CG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 91, ">=", nrow(dt_CG[dt_CG$k7_mu_rate_percentile >= 91]), nrow(dt_CG[dt_CG$k7_mu_rate_percentile < 91]))

percentile <- ecdf(dt_A$k7_mu_rate)
dt_A$k7_mu_rate_percentile <- ceiling(percentile(dt_A$k7_mu_rate)*100)
dt_A_long <- melt(dt_A, 
                      id.vars=c("k7_mu_rate", "k7_mu_rate_percentile"),
                      measure.vars=c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter"),
                      variable.name="category",
                      value.name="k7_from_N"
)
or_A_low <- odds_ratio(dt_A_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 10, "<=", nrow(dt_A[dt_A$k7_mu_rate_percentile <= 10]), nrow(dt_A[dt_A$k7_mu_rate_percentile > 10]))
or_A_high <- odds_ratio(dt_A_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 91, ">=", nrow(dt_A[dt_A$k7_mu_rate_percentile >= 91]), nrow(dt_A[dt_A$k7_mu_rate_percentile < 91]))

percentile <- ecdf(dt_C_nonCG$k7_mu_rate)
dt_C_nonCG$k7_mu_rate_percentile <- ceiling(percentile(dt_C_nonCG$k7_mu_rate)*100)
dt_C_nonCG_long <- melt(dt_C_nonCG, 
                  id.vars=c("k7_mu_rate", "k7_mu_rate_percentile"),
                  measure.vars=c("All", "CTCF binding","Enhancer - distal", "Enhancer - proximal", "Exon - CDS","Exon - UTR","Exon - other", "GERP","Miscellaneous", "Promoter"),
                  variable.name="category",
                  value.name="k7_from_N"
)
or_C_nonCG_low <- odds_ratio(dt_C_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 10, "<=", nrow(dt_C_nonCG[dt_C_nonCG$k7_mu_rate_percentile <= 10]), nrow(dt_C_nonCG[dt_C_nonCG$k7_mu_rate_percentile > 10]))
or_C_nonCG_high <- odds_ratio(dt_C_nonCG_long[,c("k7_from_N", "category", "k7_mu_rate_percentile")], 91, ">=", nrow(dt_C_nonCG[dt_C_nonCG$k7_mu_rate_percentile >= 91]), nrow(dt_C_nonCG[dt_C_nonCG$k7_mu_rate_percentile < 91]))





