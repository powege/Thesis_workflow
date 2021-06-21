rm(list = ls())
graphics.off()

library(data.table)
library(stringr)

# FUNCTION that sums across complimantary strands (ie AGC == GCT)
complement <- function(forward){
  
  complement <- stri_reverse(forward)

  # replace all bases with complement
  complement <- gsub("A", "B", complement)
  complement <- gsub("C", "D", complement)
  complement <- gsub("T", "A", complement)
  complement <- gsub("G", "C", complement)
  complement <- gsub("B", "T", complement)
  complement <- gsub("D", "G", complement)
  
return(complement)
}

# FUNCTION that takes sequence and outputs all possible sequences
all.seq.iuapc <- function(seq){
  seq <- toupper(seq)
  vec <- strsplit(seq, "")[[1]]
  vec2 <- str_replace_all(string = vec, pattern= dictio_replace)
  tmp <- expand.grid(strsplit(vec2, ""), stringsAsFactors = FALSE)
  strings <- apply(tmp, 1, paste0, collapse = "")
  return(strings)
}

# FUNCTION that gets all the 7-mers in a string
k7.window <- function(string){ mapply(function(x, y){substr(string, x, y)}, x=1:(nchar(string)-6), y=7:nchar(string)) }

### SET VARS 
dictio_replace= c("A" = "A",
                  "C" = "C",
                  "G" = "G",
                  "T" = "T",
                  "R" = "AG",        
                  "Y" = "CT",
                  "S" = "GC",
                  "W" = "AT",
                  "K" = "GT",
                  "M" = "AC",
                  "B" = "CGT",
                  "D" = "AGT",
                  "H" = "ACT",
                  "V" = "ACG",
                  "N" = "ACGT") # 

### IMPORT
htf <- fread("~/Downloads/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv")
hk7 <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/Mouse_SNV_rates/Aggarwala_Voight_2016_7mer_pmu.csv")


### FORMAT

htf <- htf[Quality == "A"]
# htf <- htf[Quality == "A" | Quality == "B"]

# length(unique(htf$`TF family`))

## TF k7
k7_list <- list()
for (j in 1:length(unique(htf$`TF family`))){
  
  sub <- htf[`TF family` == unique(htf$`TF family`)[j]]
  iuapc_strings <- sub$Consensus
  
  k7_critical_list <-  list()
  for (i in 1:length(iuapc_strings)){
    
    k7mers <- unlist(lapply(iuapc_strings[i], k7.window)) # split into 7-mers
    vec <- strsplit(iuapc_strings[i], "")[[1]] # split into vector of bases
    ind <- which(vec=="A"|vec=="C"|vec=="G"|vec=="T") # identify indicies of critical bases (ACGT)
    k7_ind <- ind[ind > 3 & ind < nchar(iuapc_strings[i])-2] -3 # calculate indicies for critical base 7-mers
    k7_critical <- k7mers[k7_ind] # subset critical base 7-mers
    k7_critical_list[[i]] <- unique(unlist(lapply(k7_critical, all.seq.iuapc))) # generate all sequences for critical 7-mers
    # print(i)
  }
  k7_list[[j]] <- unique(unlist(k7_critical_list))
  print(j)
}

htf_k7 <- sort(unique(unlist(k7_list)))
htf_k7_AC <- htf_k7[which(stri_sub(htf_k7, 4, -4) == "A" | stri_sub(htf_k7, 4, -4) == "C")]
htf_k7_GT <- htf_k7[which(stri_sub(htf_k7, 4, -4) == "G" | stri_sub(htf_k7, 4, -4) == "T")]
htf_k7 <- unique(c(htf_k7_AC, complement(htf_k7_GT)))

htf_k7 <- data.table(k7_from = htf_k7,
                     tf_status = "TF")

## k7 pmu
colnames(hk7) <- c("k7_from", "k7_to", "African_pmu", "Asian_pmu", "European_pmu", "k7_from_reverse", "k7_to_reverse")
hk7$to <- stri_sub(hk7$k7_to, 4, -4)
hk7 <- hk7[,c("k7_from", "to", "African_pmu", "Asian_pmu", "European_pmu")]

percentile <- ecdf(hk7$African_pmu)
hk7$African_pmu_percentile <-ceiling(percentile(hk7$African_pmu) * 100)


hk7 <- htf_k7[hk7, on = "k7_from"]
hk7$tf_status[is.na(hk7$tf_status)] <- "nonTF"

wilcox.test(hk7$African_pmu[hk7$tf_status == "nonTF"], hk7$African_pmu[hk7$tf_status == "TF"])
mean(hk7$African_pmu[hk7$tf_status == "nonTF"])
mean(hk7$African_pmu[hk7$tf_status == "TF"])

hk7 <- hk7[which(str_sub(hk7$k7_from, 4, -3) != "CG"),]
percentile <- ecdf(hk7$African_pmu)
hk7$African_pmu_percentile <-ceiling(percentile(hk7$African_pmu) * 100)

hist(hk7$African_pmu_percentile[hk7$tf_status == "TF"])
hist(hk7$African_pmu_percentile[hk7$tf_status == "nonTF"])

# length(unique(hk7$k7_from[hk7$tf_status == "TF"]))

A <- length(unique(hk7$k7_from[hk7$tf_status == "nonTF" &
                                 hk7$African_pmu_percentile >= 91]))
B <- length(unique(hk7$k7_from[hk7$tf_status != "nonTF" &
                                 hk7$African_pmu_percentile >= 91]))
C <- length(unique(hk7$k7_from[hk7$tf_status == "nonTF" &
                                 hk7$African_pmu_percentile < 91]))
D <- length(unique(hk7$k7_from[hk7$tf_status != "nonTF" &
                                 hk7$African_pmu_percentile < 91]))

OR <- (A/B)/(C/D)
CI <- 1.96*sqrt((1/A) + (1/B) + (1/C) + (1/D))

A <- length(unique(hk7$k7_from[hk7$tf_status == "nonTF" &
                                           hk7$African_pmu_percentile <= 20]))
B <- length(unique(hk7$k7_from[hk7$tf_status != "nonTF" &
                                           hk7$African_pmu_percentile <= 20]))
C <- length(unique(hk7$k7_from[hk7$tf_status == "nonTF" &
                                           hk7$African_pmu_percentile > 20]))
D <- length(unique(hk7$k7_from[hk7$tf_status != "nonTF" &
                                           hk7$African_pmu_percentile > 20]))

OR <- (A/B)/(C/D)
CI <- 1.96*sqrt((1/A) + (1/B) + (1/C) + (1/D))




###

# Are kmers that are more constrained in mouse than human enriched for mouse tf binding kmers that are not in human?








#####

# 
# percentile <- rep("10% most constrained", 4)
# group <- c("L", "SV", "VP", "NP")
# N <- c(length(lethal), length(subviable), length(viapheno), length(vianopheno))
# OR <- rep(NA, length(group))
# CI <- rep(NA, length(group))
# for (i in 1:length(group)){
#   A <- length(unique(dt$external_gene_name[dt$category == group[i] &
#                                              dt$Z_nonsynonymous_percentile >= 91]))
#   B <- length(unique(dt$external_gene_name[dt$category != group[i] &
#                                              dt$Z_nonsynonymous_percentile >= 91]))
#   C <- length(unique(dt$external_gene_name[dt$category == group[i] &
#                                              dt$Z_nonsynonymous_percentile < 91]))
#   D <- length(unique(dt$external_gene_name[dt$category != group[i] &
#                                              dt$Z_nonsynonymous_percentile < 91]))
#   
#   OR[i] <- (A/B)/(C/D)
#   CI[i] <- 1.96*sqrt((1/A) + (1/B) + (1/C) + (1/D))
# }
# dt_OR_zscore_top <- data.table(percentile, group, OR, CI, N)

# tf <- fread("~/Dropbox/Zhou_et_al_2017_catTFRE.csv", fill = T)
# tf$V7 <- NULL
# tf <- unlist(unlist(tf))
# tf <- tf[-which(tf=="")]
# tf <- toupper(tf)
# if (length(table(nchar(tf))) == 1){
#   tf_3.0 <- stri_sub(tf, 4, -1)
#   tf_0.3 <- stri_sub(tf, 1, -4)
#   tf_1.2 <- stri_sub(tf, 2, -3)
#   tf_2.1 <- stri_sub(tf, 3, -2)
#   tf_k7 <- c(tf_0.3, tf_3.0, tf_1.2, tf_2.1)
# }
# tf_complemnt <- stri_reverse(tf_k7)
# tf_complemnt <- gsub("A", "B", tf_complemnt)
# tf_complemnt <- gsub("C", "D", tf_complemnt)
# tf_complemnt <- gsub("T", "A", tf_complemnt)
# tf_complemnt <- gsub("G", "C", tf_complemnt)
# tf_complemnt <- gsub("B", "T", tf_complemnt)
# tf_complemnt <- gsub("D", "G", tf_complemnt)
# tf_k7 <- unique(c(tf_k7, tf_complemnt))
# 
# kmer <- fread("~/Dropbox/PhD/Data/Mu_rates/MGP_v5_allSTRAIN_kmer_PSNV_specific.csv")
# kmer <- kmer[k1_from == "A" | k1_from == "C"]
# kmer$TF <- "nonTF"
# kmer$TF[ unique(c( grep(paste(tf_k7, collapse="|"), kmer$k7_from))) ] <- "TF"
# table(kmer$TF)
# boxplot(kmer$k7_mu_rate[kmer$TF == "nonTF"], kmer$k7_mu_rate[kmer$TF == "TF"])
# hist(kmer$k7_mu_rate[kmer$TF == "nonTF"], breaks = 100)
# hist(kmer$k7_mu_rate[kmer$TF == "TF"], breaks = 100)
# mean(kmer$k7_mu_rate[kmer$TF == "nonTF"], breaks = 100)
# mean(kmer$k7_mu_rate[kmer$TF == "TF"], breaks = 100)
# 
# wilcox.test(kmer$k7_mu_rate[kmer$TF == "nonTF"], kmer$k7_mu_rate[kmer$TF == "TF"])






# k7_list <- list()
# for (j in 1:length(unique(tf$`TF family`))){
#   sub <- tf[`TF family` == unique(tf$`TF family`)[j]]
#   iuapc_strings <- sub$Consensus
#   out_list <-  list()
#   for (i in 1:length(iuapc_strings)){
#     out_list[[i]] <- all.seq.iuapc(iuapc_strings[i])
#     print(i)
#   }
#   all_strings <- unique(unlist(out_list))
#   k7_list[[j]] <- unique(unlist(lapply(all_strings, k7.window)))
#   print(j)
# }
