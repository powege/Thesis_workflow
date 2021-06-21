#!/bin/bash

# Rscript ~/Dropbox/GitHub_repos/PhD/Code/NC_constraint/Analysis/annotation_by_percentile.R \
# ~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_constraint_by_window_750_50_MGP_allMUSMUS.csv \
# ~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv \
# OER_percentile \
# 19 \
# 50 \
# ~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_annotation_by_OER_percentile_750_50_MGP_allMUSMUS.csv

# Rscript ~/Dropbox/GitHub_repos/PhD/Code/NC_constraint/Analysis/annotation_by_GW_percentile.R \
# ~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_constraint_by_window_750_50_MGP_allSTRAIN.csv \
# ~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv \
# residual_percentile \
# 19 \
# 50 \
# ~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_annotation_by_residual_percentile_750_50_MGP_allSTRAIN.csv

# Rscript ~/Dropbox/GitHub_repos/PhD/Code/NC_constraint/Analysis/annotation_by_GW_percentile.R \
# ~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_constraint_by_window_750_50_wild.csv \
# ~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv \
# OE_percentile \
# 19 \
# 50 \
# ~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_annotation_by_OE_percentile_750_50_wild.csv

Rscript ~/Dropbox/GitHub_repos/PhD/Code/NC_constraint/Analysis/annotation_by_GW_percentile.R \
~/Dropbox/PhD/Data/NC_constraint/Constraint/Mouse_constraint_by_window_750_50_wild.csv \
~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_all.csv \
residual_percentile \
19 \
50 \
~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Mouse_annotation_by_residual_percentile_750_50_wild.csv
