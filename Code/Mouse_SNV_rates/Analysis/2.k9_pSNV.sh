#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -pe shmem 10
#$ -N k9_pSNV
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

# unmethylated
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/2.k9_pSNV.R \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k9_SNV_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k5_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k9_pSNV_specific.csv.gz 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
