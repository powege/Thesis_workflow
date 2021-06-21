#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -N window_qc
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

# ### HUMAN 
# Rscript --vanilla /well/lindgren/George/Code/NC_constraint/QC/window_QC.R \
# /well/lindgren/George/Data/NC_constraint/Human_constraint_var_by_window_750_50_gnomad_MAF001_chr \
# /well/lindgren/George/Data/NC_constraint/Human_constraint_var_by_window_750_50_gnomad_MAF001_QCed.csv \
# human \
# 1 \
# 0.5

# ### MOUSE MGP MUSMUS
# Rscript --vanilla /well/lindgren/George/Code/NC_constraint/QC/window_QC.R \
# /well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_MGP_allMUSMUS_chr \
# /well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_MGP_allMUSMUS_QCed.csv \
# mouse \
# 1 \
# 0.5

### MOUSE MGP all
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/QC/window_QC.R \
/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_MGP_allSTRAIN_chr \
/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_MGP_allSTRAIN_QCed.csv \
mouse \
0.9 \
0.5

### MOUSE Wild
Rscript --vanilla /well/lindgren/George/Code/NC_constraint/QC/window_QC.R \
/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_wild_chr \
/well/lindgren/George/Data/NC_constraint/Mouse_constraint_var_by_window_750_50_wild_QCed.csv \
mouse \
0.9 \
0.5

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


