# SCRIPT that tests for enrichment of IMPC phenotypes (lethal, subviable, viable with phenotype, and viable no phenotype)
# as OR for most and least constrained genes (top and bottom 10%) (OE ratio, and nonsynonymous z-score benchmarked against spretus dNdS)
# (https://www.mousephenotype.org/impress/OntologyInfo?action=list&procID=703#31377)

rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS
odds_ratio <- function(dt_in, groups, threshold, highlow){
  
  score <- rep(colnames(dt_in)[3], length(groups))
  limit <- rep(threshold, length(groups))
  colnames(dt_in) <- c("external_gene_name", "category", "percentile")
  A <- rep(NA, length(groups))
  B <- rep(NA, length(groups))
  C <- rep(NA, length(groups))
  D <- rep(NA, length(groups))
  
  for (i in 1:length(groups)){
    if (highlow == "above"){
      A[i] <- length(unique(dt_in$external_gene_name[dt_in$category %in% groups[i] &
                                                       dt_in$percentile >= threshold]))
      B[i] <- length(unique(dt_in$external_gene_name[!dt_in$category %in% groups[i] &
                                                       dt_in$percentile >= threshold]))
      C[i] <- length(unique(dt_in$external_gene_name[dt_in$category %in% groups[i] &
                                                       dt_in$percentile < threshold]))
      D[i] <- length(unique(dt_in$external_gene_name[!dt_in$category %in% groups[i] &
                                                       dt_in$percentile < threshold]))
    }
    if (highlow == "below"){
      A[i] <- length(unique(dt_in$external_gene_name[dt_in$category %in% groups[i] &
                                                       dt_in$percentile <= threshold]))
      B[i] <- length(unique(dt_in$external_gene_name[!dt_in$category %in% groups[i] &
                                                       dt_in$percentile <= threshold]))
      C[i] <- length(unique(dt_in$external_gene_name[dt_in$category %in% groups[i] &
                                                       dt_in$percentile > threshold]))
      D[i] <- length(unique(dt_in$external_gene_name[!dt_in$category %in% groups[i] &
                                                       dt_in$percentile > threshold]))
    }
  }
  OR <- (A/B)/(C/D)
  CI95 <- 1.96*sqrt((1/A) + (1/B) + (1/C) + (1/D))
  out <- data.table(score, limit, groups, A, B, C, D, OR, CI95)
  
  return(out)
}

### SET VARS
infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv"
OR.out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_constraint_score_IMPC_viability_fertility_OR.csv"

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


## viability odds ratios

percentile_scores <- c("Z_nonsynonymous_percentile",  
                       "OE_nonsynonymous_percentile", 
                       "spretus_dNdS_percentile",
                       # "rat_dNdS_percentile",
                       "pLI_percentile",
                       "oe_lof_upper_percentile",
                       "mis_z_percentile")
out_list_high <- list()
out_list_low <- list()
for (i in 1:length(percentile_scores)){
  out_list_high[[i]] <- odds_ratio(dt_in = dt_viability[,c("mmus_external_gene_name", "viability_status", percentile_scores[i]), with = F],
                    groups = c("L", "SV", "VN", "VP"), 
                    threshold = 91, 
                    highlow = "above")
  out_list_low[[i]] <- odds_ratio(dt_in = dt_viability[,c("mmus_external_gene_name", "viability_status", percentile_scores[i]), with = F],
                     groups = c("L", "SV", "VN", "VP"), 
                     threshold = 10, 
                     highlow = "below")
}
viability_high <- do.call("rbind", out_list_high)
viability_low <- do.call("rbind", out_list_low)

## fertility odds ratios

out_list_high <- list()
out_list_low <- list()
for (i in 1:length(percentile_scores)){
  out_list_high[[i]] <- odds_ratio(dt_in = dt_fertility[,c("mmus_external_gene_name", "fertility_status", percentile_scores[i]), with = F],
                                   groups = c("FFMF", "FI", "FIMI", "MI"), 
                                   threshold = 91, 
                                   highlow = "above")
  out_list_low[[i]] <- odds_ratio(dt_in = dt_fertility[,c("mmus_external_gene_name", "fertility_status", percentile_scores[i]), with = F],
                                  groups = c("FFMF", "FI", "FIMI", "MI"), 
                                  threshold = 10, 
                                  highlow = "below")
}
fertility_high <- do.call("rbind", out_list_high)
fertility_low <- do.call("rbind", out_list_low)

viability_high$limit <- paste0(">=", viability_high$limit)
fertility_high$limit <- paste0(">=", fertility_high$limit)
viability_low$limit <- paste0("<=", viability_low$limit)
fertility_low$limit <- paste0("<=", fertility_low$limit)

output <- rbind(viability_high, viability_low, fertility_high, fertility_low)
colnames(output) <- c("score", "threshold",  "group", "A", "B", "C", "D", "OR", "CI95")

### EXPORT 
fwrite(output, OR.out.file)


#####

out2 <- output[score %in% c("OE_nonsynonymous_percentile", "spretus_dNdS_percentile") & group %in% c("L", "SV", "VP", "VN")]
out2$score[out2$score == "OE_nonsynonymous_percentile"] <- "NOER"
out2$score[out2$score == "spretus_dNdS_percentile"] <- "dNdS"
out2$OR <- round(out2$OR, 2)
out2$CI95 <- round(out2$CI95, 2)
fwrite(out2, "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Figures_tables/table_supplementary_constraint_score_IMPC_viability_OR.csv")

#####

# y <- dt_viability[OE_nonsynonymous_percentile >= 91]
# x <- y[oe_lof_upper_percentile >= 91]
# 
# 
# y <- dt_viability[OE_nonsynonymous_percentile <= 10]
# x <- y[oe_lof_upper_percentile <= 10]

# x <- dt[,c("OE_nonsynonymous", "spretus_dNdS")]
# x <- x[complete.cases(x),]
# cor.test(x$OE_nonsynonymous, x$spretus_dNdS, method = "spearman")





