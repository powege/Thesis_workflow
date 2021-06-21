#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 5
#$ -N k11_counts_null
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
# Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/4.k11_SNV_counts_by_chr_null.R \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_counts_by_chr.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_kmer_pSNV_specific.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_counts_by_chr_k11_null.csv.gz \
# /well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_counts_by_chr_k1CG_null.csv.gz 

# CG methylated
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/4.k11_SNV_counts_by_chr_null.R \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_methylated_kmer_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_methylated_k11_SNV_counts_by_chr_k11_null.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_methylated_k11_SNV_counts_by_chr_k1CG_null.csv.gz 

# CG unmethylated
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/4.k11_SNV_counts_by_chr_null.R \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_kmer_pSNV_specific.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_k11_SNV_counts_by_chr_k11_null.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_CG_unmethylated_k11_SNV_counts_by_chr_k1CG_null.csv.gz 


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


# AAAACGAAAAA,5,AAAACAAAAAA,1,15
