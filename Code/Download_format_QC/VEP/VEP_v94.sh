#!/bin/bash

while getopts i:o:s: option
do
case "${option}"
in
i) IN_FILE=${OPTARG};;
o) OUT_FILE=${OPTARG};;
s) SPECIES=${OPTARG};;
esac
done

### SCRIPT that runs Ensembl VEP v94 on formatted .vcf files (MUST BE TAB DELIMITED WITHOUT HEADER)
### ARGUMENTS: -i IN_FILE; -o OUT_FILE; -s SPECIES

# gunzip /well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allMUSMUS_snps_PASS_chr19.vcf.gz
# IN_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allMUSMUS_snps_PASS_chr19.vcf
# OUT_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/WM_Harr_etal_2016_allMUSMUS_snps_PASS_VEP_chr19.vcf
# SPECIES=mus_musculus

### load the VEP module
module load ensembl-tools/94

### set your VEP variable
VEP=/well/lindgren/George/ensembl-vep/vep

### make sure VEP version 94 is git cloned
#cd /well/lindgren/George/
#git clone https://github.com/Ensembl/ensembl-vep
#cd ensembl-vep
#git checkout release/94
#mkdir /well/lindgren/George/.vep
#perl INSTALL.pl \
#--AUTO c \
#--ASSEMBLY "GRCh38" \
#--CACHEDIR /well/lindgren/George/.vep/ \
#--SPECIES "homo_sapien" \
#--VERSION 94

# check files are in place
#ls $HOME/.vep/
#ls /well/lindgren/George/.vep/

### RUN VEP
# -i -- input file:   1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
# Flags: (https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#output)
# --cache --
# --dir_cache --
# --ccds --
# --symbol -- 
# --biotype --
# --canonical --
# --vcf -- output as vcf file
# --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" -- specify output fields
# --pick -- one annotation per variant (https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick)
$VEP --input_file "$IN_FILE" \
--dir_cache /well/lindgren/George/.vep \
--cache \
--species "$SPECIES" \
--ccds \
--symbol \
--biotype \
--canonical \
--vcf \
--fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" \
--pick \
--output_file "$OUT_FILE" \
--force_overwrite \
--no_stats \
--offline


