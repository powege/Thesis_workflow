# correlation between human and mouse pmu for nonsynsy and synonymous excluding CG>TG

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

### FUNCTIONS

# FUNCTION that returns complentary kmers (must have columns k7_from and to)
complement <- function(forward){
  
  require(stringi)
  
  complement <- forward # reverse strand
  complement$k7_from <- stri_reverse(complement$k7_from)
  complement$to <- stri_reverse(complement$to)
  
  # replace all bases with complement
  complement$k7_from <- gsub("A", "B", complement$k7_from)
  complement$k7_from <- gsub("C", "D", complement$k7_from)
  complement$k7_from <- gsub("T", "A", complement$k7_from)
  complement$k7_from <- gsub("G", "C", complement$k7_from)
  complement$k7_from <- gsub("B", "T", complement$k7_from)
  complement$k7_from <- gsub("D", "G", complement$k7_from)
  complement$to <- gsub("A", "B", complement$to)
  complement$to <- gsub("C", "D", complement$to)
  complement$to <- gsub("T", "A", complement$to)
  complement$to <- gsub("G", "C", complement$to)
  complement$to <- gsub("B", "T", complement$to)
  complement$to <- gsub("D", "G", complement$to)
  
  return(complement)
}

### SET VARS
mu.table <- "~/Dropbox/PhD/Data/Mu_rates/AA_mutation_table.csv"
k7.pSNV.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz"
hsap.k7.pSNV.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Mouse_SNV_rates/Aggarwala_Voight_2016_7mer_pmu.csv"
# out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Model_variables/WM_Harr_etal_2016_allSPECIES_transcript_k7_pSNV.csv"

### IMPORT
CT <- fread(mu.table) # mutation coding table
k7_psnv_forward <- fread(paste0("gunzip -cq ", k7.pSNV.file)) # transcript sequence
h_k7_psnv_forward <- fread(paste0("", hsap.k7.pSNV.file)) # transcript sequence

### FORMAT

# human
colnames(h_k7_psnv_forward) <- c("k7_from", "k7_to", "African_pmu", "Asian_pmu", "European_pmu", "k7_from_reverse", "k7_to_reverse")
h_k7_psnv_forward$to <- stri_sub(h_k7_psnv_forward$k7_to, 4, -4)
h_k7_psnv_forward <- h_k7_psnv_forward[,c("k7_from", "to", "African_pmu")]
h_k7_psnv <- rbind(h_k7_psnv_forward, complement(h_k7_psnv_forward)) # complement k7 pSNV

## Merge codon chnages and k7 psnv
k7_psnv <- rbind(k7_psnv_forward, complement(k7_psnv_forward)) # complement k7 pSNV
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
rm(k7_psnv_forward, k7_psnv, k7_psnv_c1, k7_psnv_c2, k7_psnv_c3, k7_psnv_c, CT, h_k7_psnv, h_k7_psnv_forward)

ct_psnv$from <- stri_sub(ct_psnv$k7_from, 4, -4)
ct_psnv_AC <- ct_psnv[from == "A" | from == "C"]
ct_psnv_GT <- complement(ct_psnv[from == "G" | from == "T"])
ct_psnv_GT$from <- stri_sub(ct_psnv_GT$k7_from, 4, -4)
ct_psnv <- rbind(ct_psnv_AC, ct_psnv_GT)
ct_psnv$mutation <- paste0(ct_psnv$from, ">", ct_psnv$to)
ct_psnv$k7_mutation <- paste0(ct_psnv$k7_from, ">", ct_psnv$to)
ct_psnv$mu_group <- ct_psnv$from
ct_psnv$mu_group[which(stri_sub(ct_psnv$k7_from, 4, -3) == "CG")] <- "CG"

# test difference between groups (all mutation types)
# test differnece between groups (A, C, CG mutations)
# test differnece between groups (A>(C,G,T), C>(A,G,T), CG>(A,G,T) mutations)


#####

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "A"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "A"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "A"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "A"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "A"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "A"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "A"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "A"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "C"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "C"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "C"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mu_group == "C"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "C"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "C"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "C"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mu_group == "C"])

###

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>C"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>C"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>C"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>C"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>C"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>C"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>C"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>C"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>T"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>T"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>T"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "A>T"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>T"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>T"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>T"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "A>T"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>A"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>A"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>A"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>A"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>A"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>A"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>A"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>A"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>G"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>G"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>G"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>G"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>G"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>G"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>G"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>G"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "non" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"])

plot(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"],
     ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"])
cor.test(ct_psnv$African_pmu[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"],
         ct_psnv$k7_mu_rate[ct_psnv$Mutation_type == "syn" & ct_psnv$mutation == "C>T" & ct_psnv$mu_group != "CG"])

x <- as.data.table(table(ct_psnv[,c("k7_from", "to", "Mutation_type")]))
table(x$N[x$Mutation_type == "non"])
y <- ct_psnv[x, on = c("k7_from", "to", "Mutation_type")]
y <- unique(y[,c("k7_from", "to", "Mutation_type", "N", "k7_mu_rate")])
y <- y[complete.cases(y),]
y_non <- y[Mutation_type == "non"]
table(y_non$N)
median(y_non$k7_mu_rate[y_non$N == 1])
median(y_non$k7_mu_rate[y_non$N == 2])
wilcox.test(y_non$k7_mu_rate[y_non$N == 1], y_non$k7_mu_rate[y_non$N == 2])



nonsense <- ct_psnv[Mutation_type == "non"]
missense <- ct_psnv[Mutation_type == "mis"]
synonymous <- ct_psnv[Mutation_type == "syn"]

non_all <- unique(nonsense[,c("k7_mutation", "mutation", "k7_mu_rate")])
mis_all <- unique(missense[,c("k7_mutation", "mutation", "k7_mu_rate")])
# mis_all <- mis_all[!k7_mutation %in% non_all$k7_mutation]
syn_all <- unique(synonymous[,c("k7_mutation", "mutation", "k7_mu_rate")])
# syn_all <- syn_all[!k7_mutation %in% c(non_all$k7_mutation)]

hist(non_all$k7_mu_rate[non_all$k7_mu_rate < 0.2])
hist(mis_all$k7_mu_rate[mis_all$k7_mu_rate < 0.2])
hist(syn_all$k7_mu_rate[syn_all$k7_mu_rate < 0.2])
median(non_all$k7_mu_rate)
median(mis_all$k7_mu_rate)
median(syn_all$k7_mu_rate)
wilcox.test(non_all$k7_mu_rate, mis_all$k7_mu_rate)$p.val
wilcox.test(non_all$k7_mu_rate, syn_all$k7_mu_rate)$p.val
wilcox.test(mis_all$k7_mu_rate, syn_all$k7_mu_rate)$p.val

non_A <- nonsense[mu_group == "A"]
non_A <- unique(non_A[,c("k7_mutation", "k7_mu_rate")])
mis_A <- missense[mu_group == "A"]
mis_A <- unique(mis_A[,c("k7_mutation", "k7_mu_rate")])
# mis_A <- mis_A[!k7_mutation %in% non_A$k7_mutation]
syn_A <- synonymous[mu_group == "A"]
syn_A <- unique(syn_A[,c("k7_mutation", "k7_mu_rate")])
# syn_A <- syn_A[!k7_mutation %in% c(non_A$k7_mutation)]

hist(non_A$k7_mu_rate)
hist(mis_A$k7_mu_rate)
hist(syn_A$k7_mu_rate)
median(non_A$k7_mu_rate)
median(mis_A$k7_mu_rate)
median(syn_A$k7_mu_rate)
wilcox.test(non_A$k7_mu_rate, mis_A$k7_mu_rate)$p.val
wilcox.test(non_A$k7_mu_rate, syn_A$k7_mu_rate)$p.val
wilcox.test(mis_A$k7_mu_rate, syn_A$k7_mu_rate)$p.val

non_C <- nonsense[from == "C"]
non_C <- unique(non_C[,c("k7_mutation", "k7_mu_rate")])
mis_C <- missense[from == "C"]
mis_C <- unique(mis_C[,c("k7_mutation", "k7_mu_rate")])
# mis_C <- mis_C[!k7_mutation %in% non_C$k7_mutation]
syn_C <- synonymous[from == "C"]
syn_C <- unique(syn_C[,c("k7_mutation", "k7_mu_rate")])
# syn_C <- syn_C[!k7_mutation %in% c(non_C$k7_mutation)]

hist(non_C$k7_mu_rate)
hist(mis_C$k7_mu_rate)
hist(syn_C$k7_mu_rate)
median(non_C$k7_mu_rate)
median(mis_C$k7_mu_rate)
median(syn_C$k7_mu_rate)
wilcox.test(non_C$k7_mu_rate, mis_C$k7_mu_rate)$p.val
wilcox.test(non_C$k7_mu_rate, syn_C$k7_mu_rate)$p.val
wilcox.test(mis_C$k7_mu_rate, syn_C$k7_mu_rate)$p.val

non_CG <- nonsense[mu_group == "CG"]
non_CG <- unique(non_CG[,c("k7_mutation", "k7_mu_rate")])
mis_CG <- missense[mu_group == "CG"]
mis_CG <- unique(mis_CG[,c("k7_mutation", "k7_mu_rate")])
# mis_CG <- mis_CG[!k7_mutation %in% non_CG$k7_mutation]
syn_CG <- synonymous[mu_group == "CG"]
syn_CG <- unique(syn_CG[,c("k7_mutation", "k7_mu_rate")])
# syn_CG <- syn_CG[!k7_mutation %in% c(non_CG$k7_mutation)]

hist(non_CG$k7_mu_rate)
hist(mis_CG$k7_mu_rate)
hist(syn_CG$k7_mu_rate)
median(non_CG$k7_mu_rate)
median(mis_CG$k7_mu_rate)
median(syn_CG$k7_mu_rate)
wilcox.test(non_CG$k7_mu_rate, mis_CG$k7_mu_rate)$p.val
wilcox.test(non_CG$k7_mu_rate, syn_CG$k7_mu_rate)$p.val
wilcox.test(mis_CG$k7_mu_rate, syn_CG$k7_mu_rate)$p.val

non_AC <- nonsense[mutation == "A>C"]
non_AC <- unique(non_AC[,c("k7_mutation", "k7_mu_rate")])
syn_AC <- synonymous[mutation == "A>C"]
syn_AC <- unique(syn_AC[,c("k7_mutation", "k7_mu_rate")])
hist(non_AC$k7_mu_rate)
hist(syn_AC$k7_mu_rate)
median(non_AC$k7_mu_rate)
median(syn_AC$k7_mu_rate)
wilcox.test(non_AC$k7_mu_rate, syn_AC$k7_mu_rate)$p.val

non_AT <- nonsense[mutation == "A>T"]
non_AT <- unique(non_AT[,c("k7_mutation", "k7_mu_rate")])
syn_AT <- synonymous[mutation == "A>T"]
syn_AT <- unique(syn_AT[,c("k7_mutation", "k7_mu_rate")])
hist(non_AT$k7_mu_rate)
hist(syn_AT$k7_mu_rate)
median(non_AT$k7_mu_rate)
median(syn_AT$k7_mu_rate)
wilcox.test(non_AT$k7_mu_rate, syn_AT$k7_mu_rate)$p.val

non_CA <- nonsense[mutation == "C>A"]
non_CA <- unique(non_CA[,c("k7_mutation", "k7_mu_rate")])
syn_CA <- synonymous[mutation == "C>A"]
syn_CA <- unique(syn_CA[,c("k7_mutation", "k7_mu_rate")])
hist(non_CA$k7_mu_rate)
hist(syn_CA$k7_mu_rate)
median(non_CA$k7_mu_rate)
median(syn_CA$k7_mu_rate)
wilcox.test(non_CA$k7_mu_rate, syn_CA$k7_mu_rate)$p.val

non_CG <- nonsense[mutation == "C>G"]
non_CG <- unique(non_CG[,c("k7_mutation", "k7_mu_rate")])
syn_CG <- synonymous[mutation == "C>G"]
syn_CG <- unique(syn_CG[,c("k7_mutation", "k7_mu_rate")])
hist(non_CG$k7_mu_rate)
hist(syn_CG$k7_mu_rate)
median(non_CG$k7_mu_rate)
median(syn_CG$k7_mu_rate)
wilcox.test(non_CG$k7_mu_rate, syn_CG$k7_mu_rate)$p.val

non_CT <- nonsense[mutation == "C>T"]
non_CT <- unique(non_CT[,c("k7_mutation", "k7_mu_rate")])
syn_CT <- synonymous[mutation == "C>T"]
syn_CT <- unique(syn_CT[,c("k7_mutation", "k7_mu_rate")])
hist(non_CT$k7_mu_rate)
hist(syn_CT$k7_mu_rate)
median(non_CT$k7_mu_rate)
median(syn_CT$k7_mu_rate)
wilcox.test(non_CT$k7_mu_rate, syn_CT$k7_mu_rate)$p.val


table(nonsense$mutation)
table(missense$mutation)
table(synonymous$mutation)
table(non_all$mutation)
table(mis_all$mutation)
table(syn_all$mutation)

AC <- ct_psnv[mutation == "A>C"]
AG <- ct_psnv[mutation == "A>G"]
AT <- ct_psnv[mutation == "A>T"]
CA <- ct_psnv[mutation == "C>A"]
CG <- ct_psnv[mutation == "C>G"]
CT <- ct_psnv[mutation == "C>T"]
CGT <- ct_psnv[mu_group == "CG" & to == "T"]

table(AC$Mutation_type)
table(AG$Mutation_type)
table(AT$Mutation_type)
table(CA$Mutation_type)
table(CG$Mutation_type)
table(CT$Mutation_type)
table(CGT$Mutation_type)




### EXPORT 
# fwrite(output, out.file)







