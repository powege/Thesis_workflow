#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-22 -tc 16
#$ -N n_mask_by_ann
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
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/Covariates/n_mask_by_annotation.R \
/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv \
/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_nMASK_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_smPOS_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_NPOS_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/gnomAD/coverage/gnomad_v2.1_fraction90_depth_less10X_chr"$SGE_TASK_ID".txt \
"$SGE_TASK_ID"

### MOUSE annotation 
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/Covariates/n_mask_by_annotation.R \
/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv \
/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_nMASK_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_smPOS_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_NPOS_Ensembl_GRC38_v94_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/MGP/Coverage/MGP_fraction90_depth_less10X_chr"$SGE_TASK_ID".txt \
"$SGE_TASK_ID"


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


