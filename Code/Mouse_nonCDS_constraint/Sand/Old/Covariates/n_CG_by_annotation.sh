#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-22 -tc 16
#$ -N n_CG_by_ann
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

# human
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/Covariates/n_CG_by_annotation.R \
"/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr""$SGE_TASK_ID"".csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_all.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_nCG_chr""$SGE_TASK_ID"".csv" \
"$SGE_TASK_ID" 

# mouse
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/Covariates/n_CG_by_annotation.R \
"/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_Ensembl_GRC38_v94_chr""$SGE_TASK_ID"".csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_nCG_chr""$SGE_TASK_ID"".csv" \
"$SGE_TASK_ID" 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
