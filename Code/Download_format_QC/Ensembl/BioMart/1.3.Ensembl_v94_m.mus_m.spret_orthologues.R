rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

### SET VARS 
out.orth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_m.spret_orthologues.csv"

### SER MARTS

# listEnsemblArchives()
# listMarts(host='http://oct2018.archive.ensembl.org')

# Rat
spret_ensembl94 <- useMart(
  host='http://oct2018.archive.ensembl.org', 
  biomart='ENSEMBL_MART_ENSEMBL', 
  dataset='mspretus_gene_ensembl'
)

# Mouse
musc_ensembl94 <- useMart(
  host='Oct2018.archive.ensembl.org', 
  biomart='ENSEMBL_MART_ENSEMBL', 
  dataset='mmusculus_gene_ensembl'
)

### PULL DATA 

# filters <- listFilters(spret_ensembl94)
# attributePages(spret_ensembl94)
# attributes <- listAttributes(spret_ensembl94, page = "homologs")
# attributes <- listAttributes(spret_ensembl94, page = "feature_page")

# Rat
spret.orth <- getBM(attributes=c(
                             'ensembl_gene_id',
                             'ensembl_transcript_id',
                             'ensembl_peptide_id',
                             'mmusculus_homolog_ensembl_gene',
                             'mmusculus_homolog_associated_gene_name',
                             'mmusculus_homolog_canonical_transcript_protein',
                             'mmusculus_homolog_ensembl_peptide',
                             'mmusculus_homolog_orthology_type',
                             'mmusculus_homolog_perc_id',
                             'mmusculus_homolog_perc_id_r1',
                             'mmusculus_homolog_goc_score',
                             'mmusculus_homolog_wga_coverage',
                             'mmusculus_homolog_dn',
                             'mmusculus_homolog_ds',
                             'mmusculus_homolog_orthology_confidence'),
                filters = c('biotype'),
                values = list("protein_coding"),
                mart = spret_ensembl94)

# Mouse
musc.orth <- getBM(attributes=c('external_gene_name',
                             'ensembl_gene_id',
                             'ensembl_transcript_id',
                             'ensembl_peptide_id',
                             'mspretus_homolog_ensembl_gene',
                             'mspretus_homolog_canonical_transcript_protein',
                             'mspretus_homolog_ensembl_peptide',
                             'mspretus_homolog_orthology_type',
                             'mspretus_homolog_perc_id',
                             'mspretus_homolog_perc_id_r1',
                             'mspretus_homolog_goc_score',
                             'mspretus_homolog_wga_coverage',
                             'mspretus_homolog_dn',
                             'mspretus_homolog_ds',
                             'mspretus_homolog_orthology_confidence'),
                filters = c('biotype'),
                values = list("protein_coding"),
                mart = musc_ensembl94)

### FORMAT

# mmusculus_homolog_canonical_transcript_protein is the canonical rat peptide ID
# remove rows with no canonical peptide ID
spret.orth <- subset(spret.orth, spret.orth$mmusculus_homolog_ensembl_gene != "")
# r.orth <- r.orth[which(r.orth$ensembl_peptide_id %in%
#                        m.orth$rnorvegicus_homolog_canonical_transcript_protein),]

musc.orth <- subset(musc.orth, musc.orth$mspretus_homolog_ensembl_gene != "")
# m.orth <- m.orth[which(m.orth$ensembl_peptide_id %in% 
#                          r.orth$mmusculus_homolog_canonical_transcript_protein),]

colnames(spret.orth) <- c(                      
                      "spretus_ensembl_gene_id",                               
                      "spretus_ensembl_transcript_id",                         
                      "spretus_ensembl_peptide_id",                            
                      "musculus_ensembl_gene_id",                
                      "musculus_external_gene_name",        
                      "mmusculus_homolog_canonical_transcript_protein",
                      "musculus_ensembl_peptide_id",             
                      "orthology_type",              
                      "Maa_match_Haa",                     
                      "Haa_match_Maa",                  
                      "goc_score",                   
                      "wga_coverage",                
                      "dn",                          
                      "ds",                          
                      "orthology_confidence")  

spret.orth <- spret.orth[c(                        
                   "spretus_ensembl_gene_id",                               
                   "spretus_ensembl_transcript_id",                         
                   "spretus_ensembl_peptide_id",                            
                   "musculus_ensembl_gene_id",                
                   "musculus_external_gene_name",        
                   "musculus_ensembl_peptide_id",             
                   "orthology_type",              
                   "Maa_match_Haa",                     
                   "Haa_match_Maa",                  
                   "goc_score",                   
                   "wga_coverage",                
                   "dn",                          
                   "ds",                          
                   "orthology_confidence")]

colnames(musc.orth) <- c("musculus_external_gene_name",                            
                      "musculus_ensembl_gene_id",                               
                      "musculus_ensembl_transcript_id",                         
                      "musculus_ensembl_peptide_id",                            
                      "spretus_ensembl_gene_id",                
                      "rnorvegicus_homolog_canonical_transcript_protein",
                      "spretus_ensembl_peptide_id",             
                      "orthology_type",              
                      "Haa_match_Maa",                     
                      "Maa_match_Haa",                  
                      "goc_score",                   
                      "wga_coverage",                
                      "dn",                          
                      "ds",                          
                      "orthology_confidence")  

musc.orth <- musc.orth[c("musculus_external_gene_name",                            
                   "musculus_ensembl_gene_id",                               
                   "musculus_ensembl_transcript_id",                         
                   "musculus_ensembl_peptide_id",                            
                   "spretus_ensembl_gene_id",                
                   "spretus_ensembl_peptide_id",             
                   "orthology_type",              
                   "Haa_match_Maa",                     
                   "Maa_match_Haa",                  
                   "goc_score",                   
                   "wga_coverage",                
                   "dn",                          
                   "ds",                          
                   "orthology_confidence")]

out <- merge(spret.orth, musc.orth, all = T)
out <- subset(out, !is.na(out$spretus_ensembl_transcript_id) & !is.na(out$musculus_ensembl_transcript_id))


### OUTPUT
fwrite(out, out.orth.file)


