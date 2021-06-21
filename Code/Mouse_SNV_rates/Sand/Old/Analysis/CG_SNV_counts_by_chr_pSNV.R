rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS
source("~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R")

### SET VARS
sm.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz"
cm.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr"
gerp.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/gerp_constrained_elements.mus_musculus.bed.gz"
ann.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
meth.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz"
snv.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_CG_snps_PASS.vcf.gz"

counts.all.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_CG_SNV_counts_by_chr.csv"
counts.meth.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_CG_methylated_SNV_counts_by_chr.csv"
counts.unmeth.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_SNV_counts_by_chr.csv"
psnv.all.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_CG_pSNV_specific.csv"
psnv.meth.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_CG_methylated_pSNV_specific.csv"
psnv.unmeth.outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_pSNV_specific.csv"


### IMPORT
sm <- fread(paste0("gunzip -cq ", sm.infile))
gerp <- fread(paste0("gunzip -cq ", gerp.infile))
cm <- list()
for (i in 1:19){
  cm[[i]] <- fread(paste0("gunzip -cq ", cm.infile, i, ".bed.gz"))
}
cm <- do.call("rbind", cm)
ann <- fread(paste0("gunzip -cq ", ann.infile))
cgs <- fread(paste0("gunzip -cq ", meth.infile))
snvs <- fread(paste0("gunzip -cq ", snv.infile))

### FORMAT

# creadt mask bed
colnames(sm) <- c("chromosome", "start", "end")
colnames(gerp) <- c("chromosome", "start", "end")
colnames(cm) <- c("chromosome", "start", "end")
ann <- ann[category %in% c("Exon - CDS", "Exon - UTR", "Exon - other", "Promoter", "Enhancer - proximal", "Enhancer - distal", "CTCF binding", "Miscellaneous", "Intron - proximal")][,c("chromosome", "start", "end")]
mask <- collapse.overlap(rbind(sm, gerp, cm, ann))

# remove masked cgs from file
ind <- bed.intersect(cgs[,c("chromosome", "start", "end")], mask) # ind is the masked CGs
ind$chromosome <- as.integer(ind$chromosome)
cgs <- cgs[!ind, on = c("chromosome", "start", "end")] # antijoinn

# subset methylated and unmethylated CGs
cg_m <- cgs[coverage > 5 & percent_methylated > 60]
cg_um <- cgs[coverage > 5 & percent_methylated < 20]
# cg_m <- cgs[-which(cgs$coverage > 5 & cgs$percent_methylated < 20)]
# cg_um <- cgs[which(cgs$coverage > 5 & cgs$percent_methylated < 20)]

# split snvs by mutation type
snv_CA = rbind( snvs[V4 == "C" & V5 == "A"], snvs[V4 == "G" & V5 == "T"])
snv_CG = rbind( snvs[V4 == "C" & V5 == "G"], snvs[V4 == "G" & V5 == "C"])
snv_CT = rbind( snvs[V4 == "C" & V5 == "T"], snvs[V4 == "G" & V5 == "A"])

# count occurences and snvs for each complementary mutation type
poliwag <- function(snv_mu, bp_to, cg_dt){
chromosome <- 1:19
from <- rep("C", 19)
to <- rep(bp_to, 19)
from_N <- rep(NA, 19)
to_N <- rep(NA, 19)
  for(i in 1:19){
    snv_sub <- snv_mu[V1 == i]
    cg_sub <- cg_dt[chromosome == i]
    from_N[i] <- nrow(cg_sub)*2
    to_N[i] <- nrow(snv_sub[V2 %in% unique(c(cg_sub$start, cg_sub$end))])
    print(i)
  }
return(data.table(chromosome, from, to, from_N, to_N))
}
counts_all <- rbind(poliwag(snv_mu = snv_CA, bp_to = "A", cg_dt = cgs),
                poliwag(snv_mu = snv_CG, bp_to = "G", cg_dt = cgs),
                poliwag(snv_mu = snv_CT, bp_to = "T", cg_dt = cgs))
counts_m <- rbind(poliwag(snv_mu = snv_CA, bp_to = "A", cg_dt = cg_m),
                poliwag(snv_mu = snv_CG, bp_to = "G", cg_dt = cg_m),
                poliwag(snv_mu = snv_CT, bp_to = "T", cg_dt = cg_m))
counts_um <- rbind(poliwag(snv_mu = snv_CA, bp_to = "A", cg_dt = cg_um),
                  poliwag(snv_mu = snv_CG, bp_to = "G", cg_dt = cg_um),
                  poliwag(snv_mu = snv_CT, bp_to = "T", cg_dt = cg_um))

psnv_all <- setDT(counts_all)[, .(from_N = sum(from_N), to_N = sum(to_N)), by = .(from, to)] # sum across chromosomes
psnv_m <- setDT(counts_m)[, .(from_N = sum(from_N), to_N = sum(to_N)), by = .(from, to)] # sum across chromosomes
psnv_um <- setDT(counts_um)[, .(from_N = sum(from_N), to_N = sum(to_N)), by = .(from, to)] # sum across chromosomes

# calculate pSNV
psnv_all$mu_rate <- psnv_all$to_N / psnv_all$from_N
psnv_m$mu_rate <- psnv_m$to_N / psnv_m$from_N
psnv_um$mu_rate <- psnv_um$to_N / psnv_um$from_N

### EXPORT
fwrite(counts_all, counts.all.outfile)
fwrite(counts_m, counts.meth.outfile)
fwrite(counts_um, counts.unmeth.outfile)
fwrite(psnv_all, psnv.all.outfile)
fwrite(psnv_m, psnv.meth.outfile)
fwrite(psnv_um, psnv.unmeth.outfile)












