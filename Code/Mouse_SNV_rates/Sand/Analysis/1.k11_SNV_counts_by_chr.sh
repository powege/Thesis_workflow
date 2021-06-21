#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 10
#$ -t 1-19 -tc 16
#$ -N k11_snv_counts
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

Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/k11_SNV_counts_by_chr.R \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr"$SGE_TASK_ID".txt.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/mmus_grcm38_v_mcar_v94_alignment_chr"$SGE_TASK_ID".long.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/mmus_grcm38_v_mpah_v94_alignment_chr"$SGE_TASK_ID".long.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allSPECIES_snps_PASS_chr"$SGE_TASK_ID".vcf.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_all_sp_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/ensembl_v94_gerp_constrained_elements.mus_musculus.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_SNV_AC_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_methylated_SNV_C_counts_by_chr.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/WM_Harr_etal_2016_allSPECIES_k11_CG_unmethylated_SNV_C_counts_by_chr.csv.gz \
"$SGE_TASK_ID" 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
