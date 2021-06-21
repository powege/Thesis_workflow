### SCRIPT that compiles data.table for analysis

rm(list=ls())
graphics.off()

library(data.table)

### SET VARS
wm.all.scores.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Constraint_scores/WM_Harr_etal_2016_allSPECIES_constraint_scores.csv"
h.scores.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Constraint_scores/gnomad.v2.1.1.lof_metrics.by_transcript.txt"
orthologues.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_m.mus_orthologues.csv"
impc.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/IMPC/Formatted/IMPC_r12_viability_fertility_pleiotropy.csv"
clinvar.infile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/Variables/ClinVar_germline_pathogenic_benign_canPC_n_SNV.csv"
outfile <- "~/Dropbox/PhD/Data/Thesis_workflow/Results/Mouse_CDS_constraint/analysis_data_table.csv"

### IMPORT
wm_scores <- fread(wm.all.scores.infile)
h_scores <- fread(h.scores.infile)
orth <- fread(orthologues.infile)
impc <- fread(impc.infile)
clinvar <- fread(clinvar.infile)

### FORMAT MOUSE SCORES
wm_scores <- wm_scores[,c("external_gene_name", 
                          "ensembl_transcript_id",
                          "cds_length",
                          "OE_nonsynonymous",
                          "Z_nonsynonymous",
                          "rat_dNdS",
                          "spretus_dNdS")]
colnames(wm_scores) <- c(c("mmus_external_gene_name", 
                           "mmus_ensembl_transcript_id",
                           "cds_length",
                           "OE_nonsynonymous",
                           "Z_nonsynonymous",
                           "rat_dNdS",
                           "spretus_dNdS"))

### FORMAT HUMAN SCORES
h_scores <- h_scores[,c("gene",
                        "transcript",
                        "pLI",
                        "oe_lof_upper",
                        "mis_z")]
colnames(h_scores) <- c("hsap_external_gene_name",
                        "hsap_ensembl_transcript_id",
                        "pLI",
                        "oe_lof_upper",
                        "mis_z")

### FORMAT ORTHOLOGUES
orth <- unique(orth[orthology_type == "ortholog_one2one"])
orth <- orth[,c("H_ensembl_transcript_id",
                "M_ensembl_transcript_id",
                "orthology_type")]
colnames(orth) <- c("hsap_ensembl_transcript_id",
                    "mmus_ensembl_transcript_id",
                    "hsap_orthology_type")

### FORMAT IMPC
names(impc)[which(names(impc) == "external_gene_name")] <- "mmus_external_gene_name"

### FORMAT ClinVar
names(clinvar)[which(names(clinvar) == "external_gene_name")] <- "hsap_external_gene_name"
names(clinvar)[which(names(clinvar) == "ensembl_transcript_id")] <- "hsap_ensembl_transcript_id"

### MERGE
dt <- merge(h_scores, orth, all = F)
dt <- merge(dt, wm_scores, all = T, by = "mmus_ensembl_transcript_id")
dt <- merge(dt, impc, all = T, by = "mmus_external_gene_name")
dt <- merge(dt, clinvar, all = T, by = c("hsap_external_gene_name", "hsap_ensembl_transcript_id"))

### PERCENTILES
percentile <- ecdf(dt$Z_nonsynonymous[!duplicated(dt$mmus_external_gene_name)])
dt$Z_nonsynonymous_percentile <- ceiling(percentile(dt$Z_nonsynonymous)*100)

percentile <- ecdf(dt$OE_nonsynonymous[!duplicated(dt$mmus_external_gene_name)])
dt$OE_nonsynonymous_percentile <- ceiling(percentile(dt$OE_nonsynonymous)*100)

percentile <- ecdf(dt$rat_dNdS[!duplicated(dt$mmus_external_gene_name)])
dt$rat_dNdS_percentile <- ceiling(percentile(dt$rat_dNdS)*100)

percentile <- ecdf(dt$spretus_dNdS[!duplicated(dt$mmus_external_gene_name)])
dt$spretus_dNdS_percentile <- ceiling(percentile(dt$spretus_dNdS)*100)

percentile <- ecdf(dt$pLI[!duplicated(dt$hsap_external_gene_name)])
dt$pLI_percentile <- ceiling(percentile(dt$pLI)*100)

percentile <- ecdf(dt$oe_lof_upper[!duplicated(dt$hsap_external_gene_name)])
dt$oe_lof_upper_percentile <- ceiling(percentile(dt$oe_lof_upper)*100)

percentile <- ecdf(dt$mis_z[!duplicated(dt$hsap_external_gene_name)])
dt$mis_z_percentile <- ceiling(percentile(dt$mis_z)*100)

### TRANSFORM RANK
cor(dt[,c("Z_nonsynonymous_percentile",
          "OE_nonsynonymous_percentile",
          "rat_dNdS_percentile",
          "spretus_dNdS_percentile",
          "pLI_percentile",
          "oe_lof_upper_percentile",
          "mis_z")],use="pairwise.complete.obs")
dt$OE_nonsynonymous_percentile <- 101 - dt$OE_nonsynonymous_percentile
dt$rat_dNdS_percentile <- 101 - dt$rat_dNdS_percentile
dt$spretus_dNdS_percentile <- 101 - dt$spretus_dNdS_percentile
dt$oe_lof_upper_percentile <- 101 - dt$oe_lof_upper_percentile

### EXPORT
fwrite(dt, outfile)

#####

# x <- dt[!is.na(dt$mmus_ensembl_transcript_id) & !is.na(dt$viability_status)]
# length(unique(x$mmus_ensembl_transcript_id))
# table(x$viability_status)



