#!/bin/bash

PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/GERP/Raw/

# wget ftp://ftp.ensembl.org/pub/release-94/bed/ensembl-compara/70_mammals.gerp_constrained_element/gerp_constrained_elements.mus_musculus.bed.gz \
# -P "$PED_ROOT"/

wget ftp://ftp.ensembl.org/pub/release-101/bed/ensembl-compara/111_mammals.gerp_constrained_element/gerp_constrained_elements.mus_musculus.bb \
-P "$PED_ROOT"/

### bigBed to bed 
module load Kent_tools/401-gompi-2019b
bigBedToBed "$PED_ROOT"gerp_constrained_elements.mus_musculus.bb "$PED_ROOT"gerp_constrained_elements.mus_musculus.bed
gzip "$PED_ROOT"gerp_constrained_elements.mus_musculus.bb
gzip "$PED_ROOT"gerp_constrained_elements.mus_musculus.bed

