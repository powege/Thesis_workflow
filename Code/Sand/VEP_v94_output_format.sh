#!/bin/bash

while getopts i:t:o: option
do
case "${option}"
in
i) IN_FILE=${OPTARG};;
t) TMP_FILE=${OPTARG};;
o) OUT_FILE=${OPTARG};;
esac
done

### SCRIPT that formats the output of VEP.sh
# subsets canonical protein coding SNVs with IMPCAT (LOW, MEDIUM, or HIGH)
### ARGUMENTS: -i IN_FILE; -t TMP_FILE; -o OUT_FILE

### set rootnames of your files
# IN_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/WM_Harr_etal_2016_allMUSMUS_snps_PASS_VEP_chr9.vcf
# TMP_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/tmp_file_WMmusmus_9
# OUT_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/WM_Harr_etal_2016_allMUSMUS_snps_PASS_VEP_canPC_IMPACT_chr9.vcf

### Save # lines 
awk '/^#/' "$IN_FILE" > "$IN_FILE".meta

### Remove # lines
sed '/^#/ d' < "$IN_FILE" > "$OUT_FILE"

### Replace '|' with '\t to split INFO column
tr '|' '\t' < "$OUT_FILE" > "$TMP_FILE" && mv "$TMP_FILE" "$OUT_FILE"

### remove "CSQ="
awk '{gsub("CSQ=", "");print}' "$OUT_FILE" > "$TMP_FILE" && mv "$TMP_FILE" "$OUT_FILE"

### Subset canonical protein-coding annotations
awk '$15 == "protein_coding" && $16 == "YES"' "$OUT_FILE" > "$TMP_FILE" && mv "$TMP_FILE" "$OUT_FILE"

### Subset IMPACT (LOW, MEDIUM, or HIGH)
awk '$12 == "LOW" ||  $12 == "MODERATE" || $12 == "HIGH"' "$OUT_FILE" > "$TMP_FILE" && mv "$TMP_FILE" "$OUT_FILE"



