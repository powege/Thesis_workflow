### TO DO
# get extended exon sequences (3 bases either side)
# run for exons and total for genes

rm(list = ls())
graphics.off()

library(data.table)
library(stringi)

### FUNCTIONS

# FUNCTION that returns complentary kmers (must have columns k7_from and to)
complement <- function(forward){
  
  complement <- forward # reverse strand
  complement$k_from <- stri_reverse(complement$k_from)
  complement$to <- stri_reverse(complement$to)
  
  # replace all bases with complement
  complement$k_from <- gsub("A", "B", complement$k_from)
  complement$k_from <- gsub("C", "D", complement$k_from)
  complement$k_from <- gsub("T", "A", complement$k_from)
  complement$k_from <- gsub("G", "C", complement$k_from)
  complement$k_from <- gsub("B", "T", complement$k_from)
  complement$k_from <- gsub("D", "G", complement$k_from)
  complement$to <- gsub("A", "B", complement$to)
  complement$to <- gsub("C", "D", complement$to)
  complement$to <- gsub("T", "A", complement$to)
  complement$to <- gsub("G", "C", complement$to)
  complement$to <- gsub("B", "T", complement$to)
  complement$to <- gsub("D", "G", complement$to)

  return(complement)
}

k9.window <- function(string){ mapply(function(x, y){substr(string, x, y)}, x=1:(nchar(string)-8), y=9:nchar(string)) }
k7.window <- function(string){ mapply(function(x, y){substr(string, x, y)}, x=1:(nchar(string)-6), y=7:nchar(string)) }
k5.window <- function(string){ mapply(function(x, y){substr(string, x, y)}, x=1:(nchar(string)-4), y=5:nchar(string)) }
k3.window <- function(string){ mapply(function(x, y){substr(string, x, y)}, x=1:(nchar(string)-2), y=3:nchar(string)) }

three.split <- function(sequence){
  vec <- strsplit(sequence, "")[[1]]
  out <- paste0(vec[c(TRUE, FALSE, FALSE)], vec[c(FALSE, TRUE, FALSE)], vec[c(FALSE, FALSE, TRUE)])
  return(out)
}

# cds_seq <- cds_seq$coding[2]
kabutops <- function(cds_seq, ct_psnv, kmer){
  
  if (kmer == 5){ 
    k5_from <- k5.window(cds_seq) 
    codons <- three.split(cds_seq) # identify codons
    codons_cut <- codons[2:(length(codons) -1)] # remove first and last codon
    Codon_from <- rep(codons_cut, times = 1, each = 3) # repeat each codon 3 times
    cds_dt <- data.table(k5_from=k5_from[2:(length(k5_from)-1)], Codon_from)
    cds_ct_psnv <- ct_psnv[cds_dt, on = c("k5_from", "Codon_from")] # merge with codon table and pSNV
    p_syn <- sum(cds_ct_psnv$k5_mu_rate[cds_ct_psnv$Mutation_type == "syn"]) # sum the probability of each mutation type
    p_mis <- sum(cds_ct_psnv$k5_mu_rate[cds_ct_psnv$Mutation_type == "mis"])
    p_non <- sum(cds_ct_psnv$k5_mu_rate[cds_ct_psnv$Mutation_type == "non"])
    output <- c(p_syn, p_mis, p_non)
    names(output) <- c("k5_p_synonymous", "k5_p_missense", "k5_p_nonsense")
    return(output)
  }
  
  if (kmer == 7){ 
    k7_from <- k7.window(cds_seq) 
  codons <- three.split(cds_seq) # identify codons
  codons_cut <- codons[2:(length(codons) -1)] # remove first and last codon
  Codon_from <- rep(codons_cut, times = 1, each = 3) # repeat each codon 3 times
  cds_dt <- data.table(k7_from, Codon_from)
  cds_ct_psnv <- ct_psnv[cds_dt, on = c("k7_from", "Codon_from")] # merge with codon table and pSNV
  p_syn <- sum(cds_ct_psnv$k7_mu_rate[cds_ct_psnv$Mutation_type == "syn"]) # sum the probability of each mutation type
  p_mis <- sum(cds_ct_psnv$k7_mu_rate[cds_ct_psnv$Mutation_type == "mis"])
  p_non <- sum(cds_ct_psnv$k7_mu_rate[cds_ct_psnv$Mutation_type == "non"])
  output <- c(p_syn, p_mis, p_non)
  names(output) <- c("k7_p_synonymous", "k7_p_missense", "k7_p_nonsense")
  return(output)
  }
}


### SET VARS
mu.table <- "~/Dropbox/PhD/Data/Mu_rates/AA_mutation_table.csv"
k7.pSNV.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/Old/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz"
k5.pSNV.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/Old/pSNV/WM_Harr_etal_2016_allSPECIES_k5_pSNV_specific.csv.gz"
cds.seq.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_canPC_seq_QCpass.csv"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/Sand/WM_Harr_etal_2016_allSPECIES_transcript_kmer_pSNV.csv"

### IMPORT
CT <- fread(mu.table) # mutation coding table
cds_seq <- fread(cds.seq.file) # transcript sequence
k7_psnv_forward <- fread(paste0("gunzip -cq ", k7.pSNV.file)) # transcript sequence
k5_psnv_forward <- fread(paste0("gunzip -cq ", k5.pSNV.file)) # transcript sequence

### FORMAT

colnames(k7_psnv_forward) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
k7_psnv <- rbind(k7_psnv_forward, complement(k7_psnv_forward)) 
colnames(k7_psnv) <- c("k7_from", "to", "k7_from_N", "k7_to_N", "k7_mu_rate")

colnames(k5_psnv_forward) <- c("k_from", "to", "k_from_N", "k_to_N", "k_mu_rate")
k5_psnv <- rbind(k5_psnv_forward, complement(k5_psnv_forward)) 
colnames(k5_psnv) <- c("k5_from", "to", "k5_from_N", "k5_to_N", "k5_mu_rate")

k7_psnv$k7_to <- paste0(substr(k7_psnv$k7_from, 1, 3), k7_psnv$to, substr(k7_psnv$k7_from, 5, 7))
k7_psnv$k5_from <- stri_sub(k7_psnv$k7_from, 2, -2)

psnv <- k7_psnv[k5_psnv, on = c("k5_from", "to")]
k7_psnv_c1 <- psnv
k7_psnv_c1$Codon_from <- substr(k7_psnv_c1$k7_from, 4, 6)
k7_psnv_c1$Codon_to <- substr(k7_psnv_c1$k7_to, 4, 6)
k7_psnv_c2 <- psnv
k7_psnv_c2$Codon_from <- substr(k7_psnv_c2$k7_from, 3, 5)
k7_psnv_c2$Codon_to <- substr(k7_psnv_c2$k7_to, 3, 5)
k7_psnv_c3 <- psnv
k7_psnv_c3$Codon_from <- substr(k7_psnv_c3$k7_from, 2, 4)
k7_psnv_c3$Codon_to <- substr(k7_psnv_c3$k7_to, 2, 4)
k7_psnv_c <- rbind(k7_psnv_c1, k7_psnv_c2, k7_psnv_c3)   
k7_psnv_c <- k7_psnv_c[order(k7_from),]
ct_psnv <- CT[k7_psnv_c, on = c("Codon_from", "Codon_to")]

## Calculate pmu
system.time({
k7_pmu <- lapply(cds_seq$coding, kabutops, ct_psnv = ct_psnv, kmer = 7) 
})
k7_pmu <- as.data.table(do.call("rbind", k7_pmu))

system.time({
  k5_pmu <- lapply(cds_seq$coding, kabutops, ct_psnv = ct_psnv, kmer = 5) 
})
k5_pmu <- as.data.table(do.call("rbind", k5_pmu))






output <- cbind(cds_seq[,c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "cds_length")], pmu)

### EXPORT 
fwrite(output, out.file)







