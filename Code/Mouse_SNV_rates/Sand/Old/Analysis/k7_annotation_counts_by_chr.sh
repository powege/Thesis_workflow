#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -pe shmem 10
#$ -t 1-22 -tc 16
#$ -N k7_ann_snv_counts
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

# mouse
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/k7_annotation_counts_by_chr.R \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr"$SGE_TASK_ID".txt.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/ensembl_v94_gerp_constrained_elements.mus_musculus.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/m.mus_GRC38_k7_counts_by_annotation_and_chr.csv.gz \
"$SGE_TASK_ID" 

# human
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_SNV_rates/Analysis/k7_annotation_counts_by_chr.R \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_seq_chr"$SGE_TASK_ID".txt.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_sm.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_N.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/ensembl_v94_gerp_constrained_elements.homo_sapiens.bed.gz \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/h.sap_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_SNV_rates/h.sap_GRC38_k7_counts_by_annotation_and_chr.csv.gz \
"$SGE_TASK_ID" 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
