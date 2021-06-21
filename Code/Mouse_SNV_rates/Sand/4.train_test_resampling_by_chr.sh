#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -pe shmem 5
#$ -t 1-100 -tc 16
#$ -N k11_resampling_by_chr
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

# k11 null
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/4.k9_train_test_resampling_by_chr.R \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr_k11_null.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_sim/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_model_fit_specific_AC_k11_null.csv \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/MAE_sim/WM_Harr_etal_2016_allSPECIES_k9_resampling_by_chr_model_fit_specific_CCG_k11_null.csv 


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

