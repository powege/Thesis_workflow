#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -pe shmem 5
#$ -t 1-100 -tc 16
#$ -N k9_resampling
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

Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/4.train_test_resampling.R \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k9_resampling_80_20_model_fit_specific_AC.csv \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k9_resampling_80_20_model_fit_specific_CCG.csv 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

