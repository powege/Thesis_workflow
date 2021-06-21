#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-19 -tc 16
#$ -N var_by_window
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

### MOUSE Wild
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_nonCDS_constraint/Variables/GW_variables_by_window.R \
"$SGE_TASK_ID" \
750 \
50 \
/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr19.txt.gz \
/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_snps_PASS_chr19.vcf.gz \
/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz \
/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz \
/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr19.bed.gz \
/Results/Mouse_SNV_rates/pSNV/WM_Harr_etal_2016_allSPECIES_k7_pSNV.csv.gz \


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


