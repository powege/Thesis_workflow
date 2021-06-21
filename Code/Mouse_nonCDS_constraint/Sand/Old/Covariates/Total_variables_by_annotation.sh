#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-22 -tc 16
#$ -N total_CG_SNV_ann_chr
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

### HUMAN annotation 
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/Covariates/Total_var_by_annotation.R \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv \
/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_controls_chr"$SGE_TASK_ID".vcf \
/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_SNV_CG_total_chr"$SGE_TASK_ID".csv \
"$SGE_TASK_ID" \
human

### MOUSE annotation 
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/Covariates/Total_var_by_annotation.R \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv \
/well/lindgren/George/Data/MGP/Variants/vcf_QCed_VEP/MGP_v5_allMUSMUS_snps_QCed_hom_chr"$SGE_TASK_ID".vcf \
/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_SNV_CG_total_chr"$SGE_TASK_ID".csv \
"$SGE_TASK_ID" \
mouse

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


