#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -pe shmem 3
#$ -t 1-22 -tc 16 
#$ -N SNV_alignnment_null
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

### LOAD MODULE
module load R/3.4.3

### SET VARS 
R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Human_mouse_mapping/Analysis/4.SNV_alignment_null.R

HUMAN_ANN_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/h.sap_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz
MOUSE_ANN_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz
ALIGNMENT_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/hsap_grch38_v_mmus_grcm38_v101_alignment_chr"$SGE_TASK_ID".bed.gz
GWAS_FILE=/well/lindgren/George/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_by_chr.csv
CHR="$SGE_TASK_ID"
FUNCTIONS_FILE=/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R
OUT_GWAS_FILE=/well/lindgren/George/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_v_mmus_RegBuild_v101_GWAS_DDS_SNV_alignment_NULL_by_chr.csv

### RUN
Rscript --vanilla $R_SCRIPT \
$HUMAN_ANN_FILE \
$MOUSE_ANN_FILE \
$ALIGNMENT_FILE \
$GWAS_FILE \
$CHR \
$FUNCTIONS_FILE \
$OUT_GWAS_FILE 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0