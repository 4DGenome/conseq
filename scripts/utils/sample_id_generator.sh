#!/bin/bash


#==================================================================================================
# Created on: 2017-04-26
# Usage: ./sample_id_generator.sh
# Author: Javier Quilez (GitHub: jaquol)
# Goal: generates unique sample identifier (SAMPLE_ID) based on the samples metadata
#==================================================================================================



#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

CONSEQ=/users/GR/mb/jquilez/projects/conseq

# variables
mode=all_samples

# paths 
url=https://zenodo.org/record/376581/files/metadata.tsv
sample_id_generator=$CONSEQ/scripts/utils/sample_id_generator.py
python=`which python`



#==================================================================================================
# COMMANDS
#==================================================================================================

# *** BEFORE PUBLICATION ***:
# uncomment `wget --no-check-certificate -q -O - $url > $spreadsheet`
# remove `cp data/metadata.tsv $spreadsheet`

# download input metadata from the above specified URL
spreadsheet=$CONSEQ/metadata/downloaded_metadata.tsv
#	wget --no-check-certificate -q -O - $url > $spreadsheet
cp metadata/metadata.tsv $spreadsheet

# generate sample IDs
echo
python $sample_id_generator $spreadsheet $mode
rm $spreadsheet
echo