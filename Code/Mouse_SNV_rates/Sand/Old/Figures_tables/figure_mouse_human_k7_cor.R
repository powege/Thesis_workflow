rm(list=ls())
graphics.off()

library(data.table)
library(ggplot2)
library(stringi)
library(MASS)

### IMPORT
mk7 <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz")
hk7 <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/Mouse_SNV_rates/Aggarwala_Voight_2016_7mer_pmu.csv")

### FORMAT
colnames(hk7) <- c("k7_from", "k7_to", "African_pmu", "Asian_pmu", "European_pmu", "k7_from_reverse", "k7_to_reverse")
hk7$to <- stri_sub(hk7$k7_to, 4, -4)
hk7 <- hk7[,c("k7_from", "to", "African_pmu", "Asian_pmu", "European_pmu")]

colnames(mk7) <- c("k7_from", "to", "k7_from_N", "k7_to_N", "mouse_pmu")
mk7 <- mk7[,c("k7_from", "to", "mouse_pmu")]

dt <- hk7[mk7, on = c("k7_from", "to")]

dt_A <- dt[which(stri_sub(dt$k7_from, 4, -4) == "A"),]
dt_C_CGT <- dt[which(stri_sub(dt$k7_from, 4, -3) == "CG" & dt$to == "T"),]
dt_C_nonCG <- dt[which(stri_sub(dt$k7_from, 4, -4) == "C" & stri_sub(dt$k7_from, 4, -3) != "CG"),]

plot(dt_A$African_pmu, dt_A$mouse_pmu)
cor.test(dt_A$African_pmu, dt_A$mouse_pmu)
plot(dt_C_nonCG$African_pmu, dt_C_nonCG$mouse_pmu)
cor.test(dt_C_nonCG$African_pmu, dt_C_nonCG$mouse_pmu)
plot(dt_C_CGT$African_pmu, dt_C_CGT$mouse_pmu)
cor.test(dt_C_CGT$African_pmu, dt_C_CGT$mouse_pmu)

dt_any <- rbind(dt_A, dt_C_nonCG)
dt_any <- setDT(dt_any)[, .(African_pmu = sum(African_pmu), mouse_pmu = sum(mouse_pmu)), by = .(k7_from)]
plot(dt_any$African_pmu, dt_any$mouse_pmu)
cor.test(dt_any$African_pmu, dt_any$mouse_pmu)




#####

# dt$mouse_pmu_rel <- (dt$mouse_pmu - mean(dt$mouse_pmu)) / mean(dt$mouse_pmu)
# dt$African_pmu_rel <- (dt$African_pmu - mean(dt$African_pmu)) / mean(dt$African_pmu)
# dt$Asian_pmu_rel <- (dt$Asian_pmu - mean(dt$Asian_pmu)) / mean(dt$Asian_pmu)
# dt$European_pmu_rel <- (dt$European_pmu - mean(dt$European_pmu)) / mean(dt$European_pmu)
# 
# 
# cor.test(dt$African_pmu_rel, dt$Asian_pmu_rel)
# cor.test(dt$African_pmu_rel, dt$European_pmu_rel)
# cor.test(dt$Asian_pmu_rel, dt$European_pmu_rel)
# plot(dt$African_pmu_rel, dt$Asian_pmu_rel)
# plot(dt$African_pmu_rel, dt$European_pmu_rel)
# plot(dt$Asian_pmu_rel, dt$European_pmu_rel)
# 
# cor.test(dt$mouse_pmu_rel, dt$African_pmu_rel)
# cor.test(dt$mouse_pmu_rel, dt$Asian_pmu_rel)
# cor.test(dt$mouse_pmu_rel, dt$European_pmu_rel)
# 
# plot(dt$mouse_pmu_rel, dt$African_pmu_rel)
# plot(dt$mouse_pmu_rel, dt$Asian_pmu_rel)
# plot(dt$mouse_pmu_rel, dt$European_pmu_rel)
# 
# plot(dt$mouse_pmu, dt$African_pmu)
# plot(dt$mouse_pmu, dt$Asian_pmu)
# plot(dt$mouse_pmu, dt$European_pmu)
# 
# # split to CG and nonCG
# dt_cgT <- dt[which(stri_sub(dt$k7_from, 4, -3) == "CG" & dt$to == "T"),]
# dt_noncgT <- dt[!dt_cgT, on = colnames(dt)]
# 
# plot(dt$mouse_pmu_rel, dt$African_pmu_rel)
# plot(dt_cgT$mouse_pmu_rel, dt_cgT$African_pmu_rel)
# plot(dt_noncgT$mouse_pmu_rel, dt_noncgT$African_pmu_rel)
# 
# cor.test(dt$mouse_pmu_rel, dt$African_pmu_rel)
# cor.test(dt_cgT$mouse_pmu_rel, dt_cgT$African_pmu_rel)
# cor.test(dt_noncgT$mouse_pmu_rel, dt_noncgT$African_pmu_rel)
# 
# mod_noncgT <- lm(mouse_pmu_rel ~ African_pmu_rel, dt_noncgT)
# summary(mod_noncgT)
# dt_noncgT$residual <- studres(mod_noncgT) 
# dt_noncgT_extremes <- dt_noncgT[residual > 5 | residual < -5]


