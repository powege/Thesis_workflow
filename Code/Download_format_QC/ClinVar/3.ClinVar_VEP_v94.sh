#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -N ClinVar_VEP
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

# Convert from .csv to .tab
tr ',' '\t' < /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps.vcf > /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps.tmp

# RUN VEP
bash /well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/VEP/VEP_v94.sh \
-i /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps.tmp \
-o /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.vcf \
-s homo_sapiens

# FORMAT VEP

# remove # lines
sed '/^#/ d' < /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.vcf > /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp

# replace '|' with '\t' to split INFO column
tr '|' '\t' < /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp > /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/tmp
mv /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/tmp /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp

# remove "CSQ="
awk '{gsub("CSQ=", "");print}' /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp > /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/tmp
mv /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/tmp /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp

# subset canonical protein-coding annotations
awk '$15 == "protein_coding" && $16 == "YES"' /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp > /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/tmp
mv /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/tmp /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp

# subset by impact
awk '$12 == "HIGH" || $12 == "MODERATE" || $12 == "LOW"' /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp > /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/tmp
mv /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/tmp /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp

# rename 
mv /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.tmp /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP_canPC_IMPACT.vcf

# gzip 
gzip /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP_canPC_IMPACT.vcf

# remove tmp files 
rm /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps.tmp
rm /well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Formatted/ClinVar_GRCh38_germline_pathogenic_benign_snps_VEP.vcf

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

## summarise consequences
#awk -F '\t' '{print $11}' M_MGP_QCed_VEP_canPC_all.txt | sort | uniq -c | sort -nr

