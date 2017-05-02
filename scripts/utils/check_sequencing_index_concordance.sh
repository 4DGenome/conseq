#!/bin/bash


#==================================================================================================
# Created on: 2017-04-21
# Usage: ./check_sequencing_index_concordance.sh
# Author: javier.quilez@crg.eu
# Goal: makes directory for analysis within the project directory
#==================================================================================================

CONSEQ=/users/GR/mb/jquilez/projects/conseq

# variables
sample_id=$1
io_metadata=$CONSEQ/scripts/utils/io_metadata.sh

# get sequencing index from the FASTQ header
# note that the index sequence may be missing from the header
ifq=data/*/raw/*/$sample_id*read1*gz
fq_header=`zcat $ifq | head -n 1 | sed s'/ /:/g'`
sequencing_index_fq=`echo $fq_header | cut -f11 -d':'`

# get sequencing index from the metadata database
sequencing_index=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a SEQUENCING_INDEX`

# Check that the sequencing index introduced as part of the metadata agrees with that found in the FASTQ reads
if [[ $sequencing_index == $sequencing_index_fq ]]; then
	echo "[OK] sequencing index is the same the metadata and FASTQ"
else
	echo "[WARN] Sequencing index added as part of the metadata does not agree with that observed in the FASTQ file"
fi
