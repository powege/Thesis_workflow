# module load R/3.4.3

rm(list = ls())
graphics.off()

library(data.table)

meth <- fread("gunzip -cq /well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz")
meth <- meth[,c("chromosome", "start", "end")]
for(i in 1:19){
  meth_sub <- meth[chromosome == i]
  vcf <- fread(paste0("gunzip -cq /well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_snps_PASS_chr", i, ".vcf.gz"))
  vcf_sub <- vcf[V2 %in% unique(c(meth_sub$start, meth_sub$end))]
  vcf_sub <- vcf_sub[,c("V1", "V2", "V4", "V5")]
  # fwrite(vcf_sub, paste0("gunzip -cq /well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_CG_snps_PASS_chr", i, ".vcf.gz"), append = T, compress = "gzip")
  fwrite(vcf_sub, "/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_CG_snps_PASS.vcf.gz", append = T, compress = "gzip")
  print(i)
  }