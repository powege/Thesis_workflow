#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -pe shmem 10
#$ -N k11_snv_counts_null
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

# all
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/3.k11_SNV_counts_by_chr_null.R \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k11_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k9_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k5_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr_k11_null.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr_k9_null.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr_k7_null.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr_k5_null.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr_k3_null.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr_k1CG_null.csv.gz

# # CG methylated
# Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/3.k11_SNV_counts_by_chr_null.R \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k9_CG_methylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_CG_methylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k5_CG_methylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_CG_methylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_CG_methylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr_k11_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr_k9_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr_k7_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr_k5_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr_k3_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr_k1CG_null.csv.gz 
# 
# # CG unmethylated
# Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/3.k11_SNV_counts_by_chr_null.R \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_obs/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_C_counts_by_chr.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k9_CG_unmethylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_CG_unmethylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k5_CG_unmethylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k3_CG_unmethylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k1CG_CG_unmethylated_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_C_counts_by_chr_k11_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_C_counts_by_chr_k9_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_C_counts_by_chr_k7_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_C_counts_by_chr_k5_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_C_counts_by_chr_k3_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/counts_sim/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_C_counts_by_chr_k1CG_null.csv.gz 


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
