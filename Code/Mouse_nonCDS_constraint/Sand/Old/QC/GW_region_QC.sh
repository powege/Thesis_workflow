#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-22 -tc 16
#$ -N region_qc
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

### HUMAN 
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/QC/region_QC.R \
/well/lindgren/George/Data/NC_constraint/Human_constraint_var_by_region_gnomad_MAF001_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/NC_constraint/Human_constraint_var_by_region_gnomad_MAF001_QCed.csv \
human \
2 \
10000 \
1 \
0.5 

### MOUSE 
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/QC/region_QC.R \
/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_region_MGP_allMUSMUS_chr"$SGE_TASK_ID".csv \
/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_region_MGP_allMUSMUS_QCed.csv \
mouse \
2 \
10000 \
1 \
0.5 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


