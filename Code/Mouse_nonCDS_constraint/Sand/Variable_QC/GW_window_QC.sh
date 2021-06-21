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

### MOUSE Wild
Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Mouse_nonCDS_constraint/Variable_QC/GW_window_QC.R \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/Variables/mmus_GRC38_constraint_var_by_winndow_750_50_WMallSP_chr \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/Variables/mmus_GRC38_constraint_var_by_window_750_50_WMallSP_QCed.csv.gz \
/well/lindgren/George/Data/Thesis_workflow/Results/Mouse_nonCDS_constraint/Variables/mmus_GRC38_constraint_var_by_window_750_50_WMallSP_removed.csv.gz \
0.9 \
0.9

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


