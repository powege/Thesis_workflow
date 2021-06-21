#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -t 1-19 -tc 16
#$ -N WM_allMUSMUS_VEP
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

# SET VARS
chromosomes=( 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X )
CHR=${chromosomes[$SGE_TASK_ID]}

# UNZIP
gunzip /well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allMUSMUS_snps_PASS_chr"$CHR".vcf.gz

# RUN VEP
bash /well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/VEP/VEP_v94.sh \
-i /well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allMUSMUS_snps_PASS_chr"$CHR".vcf \
-o /well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/WM_Harr_etal_2016_allMUSMUS_snps_PASS_VEP_chr"$CHR".vcf  \
-s mus_musculus

# FORMAT VEP
module load R/3.4.3
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/VEP/VEP_v94_output_format.R \
/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/WM_Harr_etal_2016_allMUSMUS_snps_PASS_VEP_chr"$CHR".vcf \
/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allMUSMUS_snps_PASS_VEP_canPC_IMPACT_chr"$CHR".vcf.gz

### ZIP 
gzip /well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allMUSMUS_snps_PASS_chr"$CHR".vcf
gzip /well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/WM_Harr_etal_2016_allMUSMUS_snps_PASS_VEP_chr"$CHR".vcf

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
