rm(list = ls())
graphics.off()

library(data.table)

### FUNNCTIONS
# dt_in <- dt[,c("mmus_external_gene_name","OE_nonsynonymous","n_sig_top_level_mp","n_tests_top_level_mp","hit_rate_top_level_mp")]
# minimum_tests <- 5
dugtrio <- function(dt_in, minimum_tests){
  
  colnames(dt_in) <- c("mmus_external_gene_name",
                       "score",
                       "n_sig",
                       "n_tests",
                       "hit_rate")
  dt_in <- dt_in[complete.cases(dt_in)]
  dt_in <- dt_in[n_tests >= minimum_tests] # QC minimum number of tests
  # hist(dt_in$hit_rate)
  percentile <- ecdf(dt_in$score[!duplicated(dt_in$mmus_external_gene_name)])
  dt_in$score_percentile <- 101 - (ceiling(percentile(dt_in$score)*100))
  
  # dt_mn <- aggregate(list(mean_hit_rate=dt_in$hit_rate), by=list(score_percentile=dt_in$score_percentile), FUN=mean)
  # plot(dt_mn$mean_hit_rate~dt_mn$score_percentile)
  # cor.test(dt_mn$mean_hit_rate, dt_mn$score_percentile)
  # 
  # dt_md <- aggregate(list(median_hit_rate=dt_in$hit_rate), by=list(score_percentile=dt_in$score_percentile), FUN=median)
  # plot(dt_md$median_hit_rate~dt_md$score_percentile)
  # cor.test(dt_md$median_hit_rate, dt_md$score_percentile)
  # 
  # plot(dt_in$hit_rate ~ dt_in$score_percentile)
  # mod <- glm(formula = hit_rate ~ score_percentile, data = dt_in, family = quasibinomial("logit"))
  # summary(mod)
  
  return(dt_in)
}

### SET VARS
infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv"
top.level.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_OE_IMPC_pleiotropy_top_level.csv"
procedure.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_OE_IMPC_pleiotropy_procedure.csv"
parameter.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/WM_Harr_etal_2016_allSPECIES_OE_IMPC_pleiotropy_parameter.csv"

### IMPORT
dt <- fread(infile)

### FORMAT
dt <- dt[!viability_status %in% c("L", "SV")]
dt <- dt[,c("mmus_external_gene_name",
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
            "n_sig_top_level_mp",
            "n_tests_top_level_mp",
            "n_sig_procedure",            
            "n_tests_procedure",
            "n_sig_parameter",
            "n_tests_parameter",          
            "hit_rate_top_level_mp",
            "hit_rate_procedure",
            "hit_rate_parameter")]
dt <- dt[complete.cases(dt[,c("n_sig_top_level_mp",
                              "n_tests_top_level_mp",
                              "n_sig_procedure",            
                              "n_tests_procedure",
                              "n_sig_parameter",
                              "n_tests_parameter",          
                              "hit_rate_top_level_mp",
                              "hit_rate_procedure",
                              "hit_rate_parameter")]),]

type_list <- list(c("mmus_external_gene_name","OE_nonsynonymous","n_sig_top_level_mp","n_tests_top_level_mp","hit_rate_top_level_mp"),
                  c("mmus_external_gene_name","OE_nonsynonymous","n_sig_procedure","n_tests_procedure","hit_rate_procedure"),
                  c("mmus_external_gene_name","OE_nonsynonymous","n_sig_parameter","n_tests_parameter","hit_rate_parameter"))
out_list <- list()
for (i in 1:length(type_list)){
  out_list[[i]] <- dugtrio(dt_in = dt[, c(type_list[[i]]), with = F],
                           minimum_tests = 5)
}

top_level_dt <- out_list[[1]]
procedure_dt <- out_list[[2]]
parameter_dt <- out_list[[3]]

### EXPORT 
fwrite(top_level_dt, top.level.outfile)
fwrite(procedure_dt, procedure.outfile)
fwrite(parameter_dt, parameter.outfile)


# hist(dt$hit_rate_top_level_mp)
# hist(dt$hit_rate_procedure)
# hist(dt$hit_rate_parameter)
# hist(dt$n_tests_top_level_mp)
# hist(dt$n_tests_procedure)
# hist(dt$n_tests_parameter)

#####


