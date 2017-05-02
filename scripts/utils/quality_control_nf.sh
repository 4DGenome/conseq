#!/bin/bash


#==================================================================================================
# Created on: 2017-04-07
# Usage: ./quality_control_nf.sh
# Author: Javier Quilez (GitHub: jaquol)
# Goal: perform quality control of raw reads using FastQC
#==================================================================================================

# workflow
# samples, location/type of input data, location output data/log files, cluster options
# are specified in the 'configuration variables and paths' --review before execution of the script!
# if integrate_metadata='yes':
# (1) metadata is downloaded
# (2) the FastQC output is parsed to extract metadata which is added to the database



#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

CONSEQ=/users/GR/mb/jquilez/projects/conseq

# variables
samples="gv_066_01_01_chipseq"
process=quality_control
project=jquilez
analysis=2017-04-07_analyses_manuscript
integrate_metadata="yes"
# required variables if integrate_metadat=no
release_date=
data_type=
sequencing_type=

# paths
io_metadata=/users/GR/mb/jquilez/utils/io_metadata.sh
fastqc=`which fastqc`
unzip=`which unzip`
ANALYSIS=$CONSEQ/projects/$project/analysis/$analysis	
JOB_CMD=$ANALYSIS/job_cmd
JOB_OUT=$ANALYSIS/job_out
mkdir -p $JOB_CMD
mkdir -p $JOB_OUT

# Cluster parameters
#queue=short-sl7
#memory=15G
#max_time=06:00:00
#slots=1

# download metadata and update database
if [[ $integrate_metadata == 'yes' ]]; then
	$io_metadata -m download_input_metadata
fi



#==================================================================================================
# COMMANDS
#==================================================================================================

for s in $samples; do

	# Build job: parameters
	job_name=${process}_${s}
	job_file=$JOB_CMD/$job_name.sh
	m_out=$JOB_OUT
	echo "#!/bin/bash
	#$ -N $job_name
	#$ -q $queue
	#$ -l virtual_free=$memory
	#$ -l h_rt=$max_time
	#$ -o $m_out/${job_name}_\$JOB_ID.out
	#$ -e $m_out/${job_name}_\$JOB_ID.err
	#$ -j y
	#$ -M javier.quilez@crg.eu
	#$ -m abe
	#$ -pe smp $slots" > $job_file
	sed -i 's/^\t//g' $job_file

	# FastQC
	job_cmd="$fastqc --extract $IODIR/${s}*read1.fastq.gz -o $FASTQC; rm -f $FASTQC/${s}*read1_fastqc.zip"
	echo $job_cmd >> $job_file
	if [[ $sequencing_type == "PE" ]]; then
		job_cmd="$fastqc --extract $IODIR/${s}*read2.fastq.gz -o $FASTQC; rm -f $FASTQC/${s}*read2_fastqc.zip"
		echo $job_cmd >> $job_file
	fi

	# add to metadata
	if [[ $integrate_metadata == 'yes' ]]; then
		job_cmd="$io_metadata -m quality_control_raw_reads -s $s -p $sequencing_type"
		echo $job_cmd >> $job_file
	fi

	# Submit job
	chmod a+x $job_file
	qsub < $job_file

done







IODIR=$CONSEQ/data/$data_type/raw/$release_date
FASTQC=$IODIR/fastqc
mkdir -p $FASTQC
