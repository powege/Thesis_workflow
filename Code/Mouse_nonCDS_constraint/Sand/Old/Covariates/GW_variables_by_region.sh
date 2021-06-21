#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-22 -tc 16
#$ -N var_by_region
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

echo "########################################################"
echo "Submit sample QC job"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started date: "`date`
echo "##########################################################"

module load R/3.4.3

### HUMAN 
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/Covariates/Variables_by_region.R \
human \
"$SGE_TASK_ID" \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_controls_chr"$SGE_TASK_ID".vcf \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_smPOS_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_NPOS_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/gnomAD/coverage/gnomad_v2.1_fraction90_depth_less10X_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/NC_constraint/Human_constraint_var_by_region_gnomad_MAF001_chr"$SGE_TASK_ID".csv

### MOUSE 
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/Covariates/Variables_by_region.R \
mouse \
"$SGE_TASK_ID" \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/MGP/Variants/vcf_QCed_VEP/MGP_v5_allMUSMUS_snps_QCed_hom_chr"$SGE_TASK_ID".vcf \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_smPOS_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_NPOS_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/MGP/Coverage/MGP_fraction90_depth_less10X_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_region_MGP_allMUSMUS_chr"$SGE_TASK_ID".csv

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


