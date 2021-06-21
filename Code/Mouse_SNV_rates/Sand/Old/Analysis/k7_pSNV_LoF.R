# SCRIPT that:
# classifies 7mer substitution rates by their potential to cause synonymous, missnese, and nonsesne variants.
# Tests for difference in 7mer pmu between nonsense and synonymous groups by substitution type (all substitutions, A>, C>, A>C, A>G, A>T, C>A, C>G, C>T, CG>TG). 
# Tests correlation in 7mer pmu between human and mouse for nonsense and synonymous groups by substitution type (all substitutions, A>, C>, A>C, A>G, A>T, C>A, C>G, C>T, CG>TG). 

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)
library(plyr)

### FUNCTIONS

# FUNCTION that returns complentary kmers (must have columns k7_from and to)
# forward <- k7_psnv_forward
complement <- function(forward){
  
  require(stringi)
  
  complement <- forward # reverse strand
  complement$k7_from <- stri_reverse(complement$k7_from)
  complement$k1_to <- stri_reverse(complement$k1_to)
  
  # replace all bases with complement
  complement$k7_from <- gsub("A", "B", complement$k7_from)
  complement$k7_from <- gsub("C", "D", complement$k7_from)
  complement$k7_from <- gsub("T", "A", complement$k7_from)
  complement$k7_from <- gsub("G", "C", complement$k7_from)
  complement$k7_from <- gsub("B", "T", complement$k7_from)
  complement$k7_from <- gsub("D", "G", complement$k7_from)
  complement$k1_to <- gsub("A", "B", complement$k1_to)
  complement$k1_to <- gsub("C", "D", complement$k1_to)
  complement$k1_to <- gsub("T", "A", complement$k1_to)
  complement$k1_to <- gsub("G", "C", complement$k1_to)
  complement$k1_to <- gsub("B", "T", complement$k1_to)
  complement$k1_to <- gsub("D", "G", complement$k1_to)
  
  return(complement)
}

### SET VARS
mu.table <- "~/Dropbox/PhD/Data/Mu_rates/AA_mutation_table.csv"
k7.pSNV.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV.csv.gz"
hsap.k7.pSNV.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Mouse_SNV_rates/Aggarwala_Voight_2016_7mer_pmu.csv"

ct.psnv.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k7_pSNV_codon_table.csv"
dt.counts.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k7_pSNV_codon_table_counts.csv"
dt.cor.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k7_pSNV_mmus_hsap_cor_by_mutation_effect.csv"
dt.glm.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_mutation_effect_counts~k7_pSNV.csv"
dt.wilcox.m.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k7_pSNV_mutation_effect_wilcox.csv"
dt.wilcox.h.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/hsap_k7_pSNV_mutation_effect_wilcox.csv"

### IMPORT
CT <- fread(mu.table) # mutation coding table
k7_psnv_forward <- fread(paste0("gunzip -cq ", k7.pSNV.file)) # transcript sequence
h_k7_psnv_forward <- fread(paste0("", hsap.k7.pSNV.file)) # transcript sequence

### FORMAT

# human
colnames(h_k7_psnv_forward) <- c("k7_from", "k7_to", "African_pmu", "Asian_pmu", "European_pmu", "k7_from_reverse", "k7_to_reverse")
h_k7_psnv_forward$to <- stri_sub(h_k7_psnv_forward$k7_to, 4, -4)
h_k7_psnv_forward$hsap_k7_mu_rate <- apply(h_k7_psnv_forward[,c("African_pmu", "Asian_pmu", "European_pmu")], 1, mean)
h_k7_psnv_forward <- h_k7_psnv_forward[,c("k7_from", "to", "hsap_k7_mu_rate")]
h_k7_psnv <- rbind(h_k7_psnv_forward, complement(h_k7_psnv_forward)) # complement k7 pSNV

## Merge codon chnages and k7 psnv
k7_psnv <- rbind(k7_psnv_forward[,c("k7_from", "k7_to"      "k1_from"    "k1_to"      "k7_from_N"  "k7_to_N"    "k7_mu_rate")], complement(k7_psnv_forward)) # complement k7 pSNV
k7_psnv$k7_to <- paste0(substr(k7_psnv$k7_from, 1, 3), k7_psnv$to, substr(k7_psnv$k7_from, 5, 7))
k7_psnv_c1 <- k7_psnv
k7_psnv_c1$Codon_from <- substr(k7_psnv_c1$k7_from, 4, 6)
k7_psnv_c1$Codon_to <- substr(k7_psnv_c1$k7_to, 4, 6)
k7_psnv_c2 <- k7_psnv
k7_psnv_c2$Codon_from <- substr(k7_psnv_c2$k7_from, 3, 5)
k7_psnv_c2$Codon_to <- substr(k7_psnv_c2$k7_to, 3, 5)
k7_psnv_c3 <- k7_psnv
k7_psnv_c3$Codon_from <- substr(k7_psnv_c3$k7_from, 2, 4)
k7_psnv_c3$Codon_to <- substr(k7_psnv_c3$k7_to, 2, 4)
k7_psnv_c <- rbind(k7_psnv_c1, k7_psnv_c2, k7_psnv_c3)   
k7_psnv_c <- k7_psnv_c[order(k7_from),]
ct_psnv <- CT[k7_psnv_c, on = c("Codon_from", "Codon_to")]
ct_psnv <- h_k7_psnv[ct_psnv, on = c("k7_from", "to")]
ct_psnv$from <- stri_sub(ct_psnv$k7_from, 4, -4)
ct_psnv_AC <- ct_psnv[from == "A" | from == "C"]
ct_psnv_GT <- complement(ct_psnv[from == "G" | from == "T"])
ct_psnv_GT$from <- stri_sub(ct_psnv_GT$k7_from, 4, -4)
ct_psnv <- rbind(ct_psnv_AC, ct_psnv_GT)
ct_psnv$CG_status <- "nonCG"
ct_psnv$CG_status[which(stri_sub(ct_psnv$k7_from, 4, -3) == "CG")] <- "CG"
ct_psnv$mu_group <- paste0(ct_psnv$from, ">")
# ct_psnv$mu_group[ct_psnv$CG_status == "CG" & ct_psnv$to == "T"] <- "CG>TG"
ct_psnv$mutation <- paste0(ct_psnv$from, ">", ct_psnv$to)
ct_psnv$mutation[ct_psnv$CG_status == "CG" & ct_psnv$to == "T"] <- "CG>TG"
ct_psnv$k7_mutation <- paste0(ct_psnv$k7_from, ">", ct_psnv$to)
ct_psnv$Mutation_type[ct_psnv$Mutation_type == "non"] <- "nonsense"
ct_psnv$Mutation_type[ct_psnv$Mutation_type == "mis"] <- "missense"
ct_psnv$Mutation_type[ct_psnv$Mutation_type == "syn"] <- "synonymous"
ct_psnv <- unique(ct_psnv)
# ct_psnv <- unique(ct_psnv[,c("Mutation_type", "k7_mu_rate", "mu_group",  "mutation", "k7_mutation")])
rm(complement, ct_psnv_AC, ct_psnv_GT, h_k7_psnv, h_k7_psnv_forward, hsap.k7.pSNV.file, k7_psnv,
   k7_psnv_c, k7_psnv_c1, k7_psnv_c2, k7_psnv_c3, k7_psnv_forward, k7.pSNV.file, mu.table)

### ANALYSIS

### Number of k7_substitutions by mutation type
dt_counts <- data.table(mutation=c("total", "missense", "nonsense", "synonymous"),
           n_k7_substitutions = c(length(unique(ct_psnv$k7_mutation)),
                                  nrow(unique(ct_psnv[Mutation_type == "missense"][,c("k7_mutation")])),
                                  nrow(unique(ct_psnv[Mutation_type == "synonymous"][,c("k7_mutation")])),
                                  nrow(unique(ct_psnv[Mutation_type == "nonsense"][,c("k7_mutation")]))))


### Test relationship between the number of potential missense mutation and k7 substitution rate. 
out_list <- list()
for(i in 1:length(unique(ct_psnv$Mutation_type))){
  x <- as.data.table(table(ct_psnv[Mutation_type == unique(ct_psnv$Mutation_type)[i]][,c("k7_mutation")]))
  colnames(x) <- c("k7_mutation", "N")
  x <- ct_psnv[,c("k7_mutation", "k7_mu_rate")][x, on = "k7_mutation"]
  x <- unique(x)
  # cor <- cor.test(x$k7_mu_rate, x$N, method = "spearman")
  mod <- glm(x$N ~ x$k7_mu_rate, family = "poisson")
  out_list[[i]] <- c(Mutation_type=unique(ct_psnv$Mutation_type)[i], 
                     mu_rate_N=nrow(x),
                     potential_N=max(x$N),
                     coef=coef(summary(mod))[2,1],
                     p_val=coef(summary(mod))[2,4])
}
dt_glm <- as.data.table(do.call("rbind", out_list))

### Test correlation between huaman and mouse substitution rates by mutation type
tmp2 <- ct_psnv
tmp2$mutation <- "All"
dt <- rbind(ct_psnv, tmp2)
# dt <- unique(dt[,c("k7_mutation", "mutation", "Mutation_type", "hsap_k7_mu_rate", "k7_mu_rate")])

out_list <- list()
for (i in 1:length(unique(dt$mutation))){
  sub <- dt[mutation == unique(dt$mutation)[i]]
  out_list2 <- list()
  for (j in 1:length(unique(sub$Mutation_type))){
    sub2 <- sub[Mutation_type == unique(sub$Mutation_type)[j]]
    cor <- cor.test(sub2$hsap_k7_mu_rate, sub2$k7_mu_rate, method = "spearman")
    out_list2[[j]] <- c(Mutation_group=unique(dt$mutation)[i],
                        Mutation_type=unique(sub$Mutation_type)[j],
                        N=nrow(sub2), 
                        rho=cor$estimate, 
                        p.val=cor$p.value)
  }
  out_list[[i]] <- do.call("rbind", out_list2)
}
dt_cor <- as.data.table(do.call("rbind", out_list))

### Test differneces in distributions of substitution probabilities by mutation type (Wilcox test)

## MOUSE
tmp2 <- ct_psnv
tmp2$mutation <- "All"
dt <- rbind(ct_psnv, tmp2)
# dt <- unique(dt[,c("k7_mutation", "mutation", "Mutation_type", "hsap_k7_mu_rate", "k7_mu_rate")])

# synonymous vs nonsense
syn_non_list <- list()
for (i in 1:length(unique(dt$mutation))){
  sub <- dt[mutation == unique(dt$mutation)[i]]
  if(unique(sub$mutation) != "A>G"){
    w_syn_non <- wilcox.test(sub$k7_mu_rate[sub$Mutation_type == "synonymous"],
                             sub$k7_mu_rate[sub$Mutation_type == "nonsense"])
    syn_non_list[[i]] <- c(Mutation_group=unique(dt$mutation)[i],
                           synonymous_N=nrow(sub[Mutation_type == "synonymous"]),
                           synonymous_median=median(sub$k7_mu_rate[sub$Mutation_type == "synonymous"]),
                           nonsense_N=nrow(sub[Mutation_type == "nonsense"]),
                           nonsense_median=median(sub$k7_mu_rate[sub$Mutation_type == "nonsense"]),
                           statistic=w_syn_non$statistic,
                           p_val=w_syn_non$p.value)
    
  }
}
dt_syn_non <- as.data.table(do.call("rbind", syn_non_list))

# synonymous vs missense
syn_mis_list <- list()
for (i in 1:length(unique(dt$mutation))){
  sub <- dt[mutation == unique(dt$mutation)[i]]
  w_syn_mis <- wilcox.test(sub$k7_mu_rate[sub$Mutation_type == "synonymous"],
                           sub$k7_mu_rate[sub$Mutation_type == "missense"])
  syn_mis_list[[i]] <- c(Mutation_group=unique(dt$mutation)[i],
                         synonymous_N=nrow(sub[Mutation_type == "synonymous"]),
                         synonymous_median=median(sub$k7_mu_rate[sub$Mutation_type == "synonymous"]),
                         missense_N=nrow(sub[Mutation_type == "missense"]),
                         missense_median=median(sub$k7_mu_rate[sub$Mutation_type == "missense"]),
                         statistic=w_syn_mis$statistic,
                         p_val=w_syn_mis$p.value)
  
}
dt_syn_mis <- as.data.table(do.call("rbind", syn_mis_list))

# nonsense vs missense
non_mis_list <- list()
for (i in 1:length(unique(dt$mutation))){
  sub <- dt[mutation == unique(dt$mutation)[i]]
  if(unique(sub$mutation) != "A>G"){
    w_non_mis <- wilcox.test(sub$k7_mu_rate[sub$Mutation_type == "nonsense"],
                             sub$k7_mu_rate[sub$Mutation_type == "missense"])
    non_mis_list[[i]] <- c(Mutation_group=unique(dt$mutation)[i],
                           nonsense_N=nrow(sub[Mutation_type == "nonsense"]),
                           nonsense_median=median(sub$k7_mu_rate[sub$Mutation_type == "nonsense"]),
                           missense_N=nrow(sub[Mutation_type == "missense"]),
                           missense_median=median(sub$k7_mu_rate[sub$Mutation_type == "missense"]),
                           statistic=w_non_mis$statistic,
                           p_val=w_non_mis$p.value)
    
  }
}
dt_non_mis <- as.data.table(do.call("rbind", non_mis_list))

dt_wilcox_m <- rbind.fill(dt_non_mis, dt_syn_non, dt_syn_mis)
dt_wilcox_m <- dt_wilcox_m[,c("Mutation_group", 
                          "nonsense_N", "nonsense_median", 
                          "missense_N", "missense_median",
                          "synonymous_N", "synonymous_median",
                          "statistic.W", "p_val")]

## HUMAN
tmp2 <- ct_psnv
tmp2$mutation <- "All"
dt <- rbind(ct_psnv, tmp2)
# dt <- unique(dt[,c("k7_mutation", "mutation", "Mutation_type", "hsap_k7_mu_rate", "k7_mu_rate")])

# synonymous vs nonsense
syn_non_list <- list()
for (i in 1:length(unique(dt$mutation))){
  sub <- dt[mutation == unique(dt$mutation)[i]]
  if(unique(sub$mutation) != "A>G"){
    w_syn_non <- wilcox.test(sub$hsap_k7_mu_rate[sub$Mutation_type == "synonymous"],
                             sub$hsap_k7_mu_rate[sub$Mutation_type == "nonsense"])
    syn_non_list[[i]] <- c(Mutation_group=unique(dt$mutation)[i],
                           synonymous_N=nrow(sub[Mutation_type == "synonymous"]),
                           synonymous_median=median(sub$hsap_k7_mu_rate[sub$Mutation_type == "synonymous"]),
                           nonsense_N=nrow(sub[Mutation_type == "nonsense"]),
                           nonsense_median=median(sub$hsap_k7_mu_rate[sub$Mutation_type == "nonsense"]),
                           statistic=w_syn_non$statistic,
                           p_val=w_syn_non$p.value)
    
  }
}
dt_syn_non <- as.data.table(do.call("rbind", syn_non_list))

# synonymous vs missense
syn_mis_list <- list()
for (i in 1:length(unique(dt$mutation))){
  sub <- dt[mutation == unique(dt$mutation)[i]]
  w_syn_mis <- wilcox.test(sub$hsap_k7_mu_rate[sub$Mutation_type == "synonymous"],
                           sub$hsap_k7_mu_rate[sub$Mutation_type == "missense"])
  syn_mis_list[[i]] <- c(Mutation_group=unique(dt$mutation)[i],
                         synonymous_N=nrow(sub[Mutation_type == "synonymous"]),
                         synonymous_median=median(sub$hsap_k7_mu_rate[sub$Mutation_type == "synonymous"]),
                         missense_N=nrow(sub[Mutation_type == "missense"]),
                         missense_median=median(sub$hsap_k7_mu_rate[sub$Mutation_type == "missense"]),
                         statistic=w_syn_mis$statistic,
                         p_val=w_syn_mis$p.value)
  
}
dt_syn_mis <- as.data.table(do.call("rbind", syn_mis_list))

# nonsense vs missense
non_mis_list <- list()
for (i in 1:length(unique(dt$mutation))){
  sub <- dt[mutation == unique(dt$mutation)[i]]
  if(unique(sub$mutation) != "A>G"){
    w_non_mis <- wilcox.test(sub$hsap_k7_mu_rate[sub$Mutation_type == "nonsense"],
                             sub$hsap_k7_mu_rate[sub$Mutation_type == "missense"])
    non_mis_list[[i]] <- c(Mutation_group=unique(dt$mutation)[i],
                           nonsense_N=nrow(sub[Mutation_type == "nonsense"]),
                           nonsense_median=median(sub$hsap_k7_mu_rate[sub$Mutation_type == "nonsense"]),
                           missense_N=nrow(sub[Mutation_type == "missense"]),
                           missense_median=median(sub$hsap_k7_mu_rate[sub$Mutation_type == "missense"]),
                           statistic=w_non_mis$statistic,
                           p_val=w_non_mis$p.value)
    
  }
}
dt_non_mis <- as.data.table(do.call("rbind", non_mis_list))

dt_wilcox_h <- rbind.fill(dt_non_mis, dt_syn_non, dt_syn_mis)
dt_wilcox_h <- dt_wilcox_h[,c("Mutation_group", 
                          "nonsense_N", "nonsense_median", 
                          "missense_N", "missense_median",
                          "synonymous_N", "synonymous_median",
                          "statistic.W", "p_val")]

### EXPORT

fwrite(ct_psnv, ct.psnv.outfile)
fwrite(dt_counts, dt.counts.outfile)
fwrite(dt_cor, dt.cor.outfile)
fwrite(dt_glm, dt.glm.outfile)
fwrite(dt_wilcox_m, dt.wilcox.m.outfile)
fwrite(dt_wilcox_h, dt.wilcox.h.outfile)

####


# p_box <- ggplot(ct_psnv, mapping = aes(x=mutation, y=k7_mu_rate)) +
#   # geom_violin(width=1) +
#   geom_boxplot(width=0.2, outlier.shape = NA, fill = "white") +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   # ggtitle("") +
#   ylim(0, 0.1) +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   )
# p_box
# 
# p_box_C <- ggplot(ct_psnv[mu_group == "C>"], mapping = aes(x=mutation, y=k7_mu_rate)) +
#   geom_violin(width=1.2) +
#   geom_boxplot(width=0.15, outlier.shape = NA, fill = "white") +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   ylim(0, 0.2) +
#   # ggtitle("") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# p_box_C
# 
# p_box_A <- ggplot(ct_psnv[mu_group == "A>"], mapping = aes(x=mutation, y=k7_mu_rate)) +
#   geom_violin(width=1.2) +
#   geom_boxplot(width=0.15, outlier.shape = NA, fill = "white") +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   # ggtitle("") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# p_box_A
# 
# p_box_CGTG <- ggplot(ct_psnv[mu_group == "CG>TG"], mapping = aes(x=mutation, y=k7_mu_rate)) +
#   # geom_violin(width=1.2) +
#   geom_boxplot(width=0.15, outlier.shape = NA, fill = "white") +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   # ggtitle("") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# p_box_CGTG
# 
# tmp <- ct_psnv[Mutation_type %in% c("nonsense", "synonymous", "missense") & mu_group %in% c("A>", "C>")]
# tmp$mu_group <- "(A|C)>"
# dt_plot <- rbind(tmp, ct_psnv[Mutation_type %in% c("nonsense", "synonymous", "missense") & mu_group %in% c("A>", "C>")])
# 
# # p_box <- ggplot(dt_plot, mapping = aes(x=mutation, y=k7_mu_rate, fill=Mutation_type)) +
# #   geom_boxplot(width=0.6, outlier.shape = NA, position=position_dodge(0.8)) +
# #   scale_fill_manual(values=c("#0072B2", "#D55E00")) +
# #   # geom_hline(yintercept=50, linetype="dashed", color = "black", size = 0.9) +
# #   ylab("7-mer substitution rate") +
# #   xlab("") +
# #   # ggtitle("(A) Distributions") +
# #   ylim(0, 0.1) +
# #   # scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.15)) +
# #   coord_flip() +
# #   theme_classic() +
# #   theme(axis.text = element_text(size = 14),
# #         axis.title = element_text(size = 14),
# #         plot.title =  element_text(size = 20, face = "bold"),
# #         # legend.position = "top",
# #         legend.title = element_blank(),
# #         legend.text = element_text(size = 14),
# #         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
# #         panel.border = element_rect(colour = "black", fill=NA, size=1)
# #   )
# # p_box
# 
# p_box <- ggplot(dt_plot, mapping = aes(x=mu_group, y=k7_mu_rate, fill=Mutation_type)) +
#   geom_boxplot(width=0.6, outlier.shape = NA, position=position_dodge(0.8)) +
#   # scale_fill_manual(values=c("#0072B2", "#D55E00")) +
#   # geom_hline(yintercept=50, linetype="dashed", color = "black", size = 0.9) +
#   ylab("7-mer substitution rate") +
#   xlab("") +
#   # ggtitle("(A) Distributions") +
#   ylim(0, 0.15) +
#   # scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.15)) +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         # legend.position = "top",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 14),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# p_box
# 
# dt_counts <- as.data.table(table(ct_psnv[,c("mutation", "Mutation_type")]))
# # dt_counts <- dt_counts[Mutation_type != "missense"]
# dt_counts <- dt_counts[, total:=sum(N), by=mutation]
# dt_counts$fraction <- dt_counts$N / dt_counts$total
# # dt_counts <- dt_counts[Mutation_type != "missense"]
# 
# ggplot(data=dt_counts, aes(x=mutation, y=N, fill=Mutation_type)) +
#   geom_bar(stat="identity", colour="black") +
#   coord_flip() +
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         plot.title =  element_text(size = 20, face = "bold"),
#         # legend.position = "top",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 14),
#         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1)
#   ) 
# 
# # dt_counts <- as.data.table(table(ct_psnv[,c("mutation", "Mutation_type")]))
# # dt_counts <- dt_counts[Mutation_type != "missense"]
# # dt_counts <- dt_counts[, total:=sum(N), by=Mutation_type]
# # dt_counts$fraction <- dt_counts$N / dt_counts$total
# # # dt_counts <- dt_counts[Mutation_type != "missense"]
# # 
# # ggplot(data=dt_counts, aes(x=Mutation_type, y=fraction, fill=mutation)) +
# #   geom_bar(stat="identity", colour="black") +
# #   coord_flip() +
# #   theme_classic() +
# #   theme(axis.text = element_text(size = 14),
# #         axis.title = element_text(size = 14),
# #         plot.title =  element_text(size = 20, face = "bold"),
# #         # legend.position = "top",
# #         legend.title = element_blank(),
# #         legend.text = element_text(size = 14),
# #         plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
# #         panel.border = element_rect(colour = "black", fill=NA, size=1)
# #   ) 
# 
# #####
# 
# 
# 
# #####
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "A"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "A"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "A"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "A"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "A"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "A"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "A"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "A"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "C"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "C"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "C"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "C"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "C"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "C"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "C"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "C"])
# 
# ###
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>C"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>C"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>C"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>C"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>C"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>C"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>C"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>C"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>T"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>T"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>T"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>T"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>T"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>T"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>T"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>T"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>A"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>A"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>A"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>A"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>A"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>A"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>A"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>A"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>G"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>G"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>G"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>G"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>G"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>G"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>G"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>G"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"])
# 
# plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"],
#      ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"])
# cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"],
#          ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"])
# 
# x <- as.data.table(table(ct_psnv[,c("k7_from", "to", "Mutation_type")]))
# table(x$N[x$Mutation_type == "non"])
# y <- ct_psnv[x, on = c("k7_from", "to", "Mutation_type")]
# y <- unique(y[,c("k7_from", "to", "Mutation_type", "N", "k7_mu_rate")])
# y <- y[complete.cases(y),]
# y_non <- y[Mutation_type == "non"]
# table(y_non$N)
# median(y_non$k7_mu_rate[y_non$N == 1])
# median(y_non$k7_mu_rate[y_non$N == 2])
# wilcox.test(y_non$k7_mu_rate[y_non$N == 1], y_non$k7_mu_rate[y_non$N == 2])
# 
# 
# 
# nonsense <- ct_psnv[Mutation_type == "non"]
# missense <- ct_psnv[Mutation_type == "mis"]
# synonymous <- ct_psnv[Mutation_type == "syn"]
# 
# non_all <- unique(nonsense[,c("k7_mutation", "mutation", "k7_mu_rate")])
# mis_all <- unique(missense[,c("k7_mutation", "mutation", "k7_mu_rate")])
# # mis_all <- mis_all[!k7_mutation %in% non_all$k7_mutation]
# syn_all <- unique(synonymous[,c("k7_mutation", "mutation", "k7_mu_rate")])
# # syn_all <- syn_all[!k7_mutation %in% c(non_all$k7_mutation)]
# 
# hist(non_all$k7_mu_rate[non_all$k7_mu_rate < 0.2])
# hist(mis_all$k7_mu_rate[mis_all$k7_mu_rate < 0.2])
# hist(syn_all$k7_mu_rate[syn_all$k7_mu_rate < 0.2])
# median(non_all$k7_mu_rate)
# median(mis_all$k7_mu_rate)
# median(syn_all$k7_mu_rate)
# wilcox.test(non_all$k7_mu_rate, mis_all$k7_mu_rate)$p.val
# wilcox.test(non_all$k7_mu_rate, syn_all$k7_mu_rate)$p.val
# wilcox.test(mis_all$k7_mu_rate, syn_all$k7_mu_rate)$p.val
# 
# non_A <- nonsense[mu_group == "A"]
# non_A <- unique(non_A[,c("k7_mutation", "k7_mu_rate")])
# mis_A <- missense[mu_group == "A"]
# mis_A <- unique(mis_A[,c("k7_mutation", "k7_mu_rate")])
# # mis_A <- mis_A[!k7_mutation %in% non_A$k7_mutation]
# syn_A <- synonymous[mu_group == "A"]
# syn_A <- unique(syn_A[,c("k7_mutation", "k7_mu_rate")])
# # syn_A <- syn_A[!k7_mutation %in% c(non_A$k7_mutation)]
# 
# hist(non_A$k7_mu_rate)
# hist(mis_A$k7_mu_rate)
# hist(syn_A$k7_mu_rate)
# median(non_A$k7_mu_rate)
# median(mis_A$k7_mu_rate)
# median(syn_A$k7_mu_rate)
# wilcox.test(non_A$k7_mu_rate, mis_A$k7_mu_rate)$p.val
# wilcox.test(non_A$k7_mu_rate, syn_A$k7_mu_rate)$p.val
# wilcox.test(mis_A$k7_mu_rate, syn_A$k7_mu_rate)$p.val
# 
# non_C <- nonsense[from == "C"]
# non_C <- unique(non_C[,c("k7_mutation", "k7_mu_rate")])
# mis_C <- missense[from == "C"]
# mis_C <- unique(mis_C[,c("k7_mutation", "k7_mu_rate")])
# # mis_C <- mis_C[!k7_mutation %in% non_C$k7_mutation]
# syn_C <- synonymous[from == "C"]
# syn_C <- unique(syn_C[,c("k7_mutation", "k7_mu_rate")])
# # syn_C <- syn_C[!k7_mutation %in% c(non_C$k7_mutation)]
# 
# hist(non_C$k7_mu_rate)
# hist(mis_C$k7_mu_rate)
# hist(syn_C$k7_mu_rate)
# median(non_C$k7_mu_rate)
# median(mis_C$k7_mu_rate)
# median(syn_C$k7_mu_rate)
# wilcox.test(non_C$k7_mu_rate, mis_C$k7_mu_rate)$p.val
# wilcox.test(non_C$k7_mu_rate, syn_C$k7_mu_rate)$p.val
# wilcox.test(mis_C$k7_mu_rate, syn_C$k7_mu_rate)$p.val
# 
# non_CG <- nonsense[mu_group == "CG"]
# non_CG <- unique(non_CG[,c("k7_mutation", "k7_mu_rate")])
# mis_CG <- missense[mu_group == "CG"]
# mis_CG <- unique(mis_CG[,c("k7_mutation", "k7_mu_rate")])
# # mis_CG <- mis_CG[!k7_mutation %in% non_CG$k7_mutation]
# syn_CG <- synonymous[mu_group == "CG"]
# syn_CG <- unique(syn_CG[,c("k7_mutation", "k7_mu_rate")])
# # syn_CG <- syn_CG[!k7_mutation %in% c(non_CG$k7_mutation)]
# 
# hist(non_CG$k7_mu_rate)
# hist(mis_CG$k7_mu_rate)
# hist(syn_CG$k7_mu_rate)
# median(non_CG$k7_mu_rate)
# median(mis_CG$k7_mu_rate)
# median(syn_CG$k7_mu_rate)
# wilcox.test(non_CG$k7_mu_rate, mis_CG$k7_mu_rate)$p.val
# wilcox.test(non_CG$k7_mu_rate, syn_CG$k7_mu_rate)$p.val
# wilcox.test(mis_CG$k7_mu_rate, syn_CG$k7_mu_rate)$p.val
# 
# non_AC <- nonsense[mutation == "A>C"]
# non_AC <- unique(non_AC[,c("k7_mutation", "k7_mu_rate")])
# syn_AC <- synonymous[mutation == "A>C"]
# syn_AC <- unique(syn_AC[,c("k7_mutation", "k7_mu_rate")])
# hist(non_AC$k7_mu_rate)
# hist(syn_AC$k7_mu_rate)
# median(non_AC$k7_mu_rate)
# median(syn_AC$k7_mu_rate)
# wilcox.test(non_AC$k7_mu_rate, syn_AC$k7_mu_rate)$p.val
# 
# non_AT <- nonsense[mutation == "A>T"]
# non_AT <- unique(non_AT[,c("k7_mutation", "k7_mu_rate")])
# syn_AT <- synonymous[mutation == "A>T"]
# syn_AT <- unique(syn_AT[,c("k7_mutation", "k7_mu_rate")])
# hist(non_AT$k7_mu_rate)
# hist(syn_AT$k7_mu_rate)
# median(non_AT$k7_mu_rate)
# median(syn_AT$k7_mu_rate)
# wilcox.test(non_AT$k7_mu_rate, syn_AT$k7_mu_rate)$p.val
# 
# non_CA <- nonsense[mutation == "C>A"]
# non_CA <- unique(non_CA[,c("k7_mutation", "k7_mu_rate")])
# syn_CA <- synonymous[mutation == "C>A"]
# syn_CA <- unique(syn_CA[,c("k7_mutation", "k7_mu_rate")])
# hist(non_CA$k7_mu_rate)
# hist(syn_CA$k7_mu_rate)
# median(non_CA$k7_mu_rate)
# median(syn_CA$k7_mu_rate)
# wilcox.test(non_CA$k7_mu_rate, syn_CA$k7_mu_rate)$p.val
# 
# non_CG <- nonsense[mutation == "C>G"]
# non_CG <- unique(non_CG[,c("k7_mutation", "k7_mu_rate")])
# syn_CG <- synonymous[mutation == "C>G"]
# syn_CG <- unique(syn_CG[,c("k7_mutation", "k7_mu_rate")])
# hist(non_CG$k7_mu_rate)
# hist(syn_CG$k7_mu_rate)
# median(non_CG$k7_mu_rate)
# median(syn_CG$k7_mu_rate)
# wilcox.test(non_CG$k7_mu_rate, syn_CG$k7_mu_rate)$p.val
# 
# non_CT <- nonsense[mutation == "C>T"]
# non_CT <- unique(non_CT[,c("k7_mutation", "k7_mu_rate")])
# syn_CT <- synonymous[mutation == "C>T"]
# syn_CT <- unique(syn_CT[,c("k7_mutation", "k7_mu_rate")])
# hist(non_CT$k7_mu_rate)
# hist(syn_CT$k7_mu_rate)
# median(non_CT$k7_mu_rate)
# median(syn_CT$k7_mu_rate)
# wilcox.test(non_CT$k7_mu_rate, syn_CT$k7_mu_rate)$p.val
# 
# 
# table(nonsense$mutation)
# table(missense$mutation)
# table(synonymous$mutation)
# table(non_all$mutation)
# table(mis_all$mutation)
# table(syn_all$mutation)
# 
# AC <- ct_psnv[mutation == "A>C"]
# AG <- ct_psnv[mutation == "A>G"]
# AT <- ct_psnv[mutation == "A>T"]
# CA <- ct_psnv[mutation == "C>A"]
# CG <- ct_psnv[mutation == "C>G"]
# CT <- ct_psnv[mutation == "C>T"]
# CGT <- ct_psnv[mu_group == "CG" & to == "T"]
# 
# table(AC$Mutation_type)
# table(AG$Mutation_type)
# table(AT$Mutation_type)
# table(CA$Mutation_type)
# table(CG$Mutation_type)
# table(CT$Mutation_type)
# table(CGT$Mutation_type)




### EXPORT 
# fwrite(output, out.file)







