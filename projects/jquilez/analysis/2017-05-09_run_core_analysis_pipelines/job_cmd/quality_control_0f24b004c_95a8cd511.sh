#!/bin/bash
#$ -N quality_control_0f24b004c_95a8cd511
#$ -q short-sl7
#$ -l virtual_free=15G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/projects/conseq/projects/jquilez/analysis/2017-05-09_run_core_analysis_pipelines/job_out/quality_control_0f24b004c_95a8cd511_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/conseq/projects/jquilez/analysis/2017-05-09_run_core_analysis_pipelines/job_out/quality_control_0f24b004c_95a8cd511_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/conseq/data/hic/raw/2015-08-01/0f24b004c_95a8cd511_read1.fastq.gz -o /users/GR/mb/jquilez/projects/conseq/data/hic/raw/2015-08-01/fastqc; rm -f /users/GR/mb/jquilez/projects/conseq/data/hic/raw/2015-08-01/fastqc/0f24b004c_95a8cd511*read1_fastqc.zip
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/conseq/data/hic/raw/2015-08-01/0f24b004c_95a8cd511_read2.fastq.gz -o /users/GR/mb/jquilez/projects/conseq/data/hic/raw/2015-08-01/fastqc; rm -f /users/GR/mb/jquilez/projects/conseq/data/hic/raw/2015-08-01/fastqc/0f24b004c_95a8cd511*read2_fastqc.zip
/users/GR/mb/jquilez/projects/conseq/scripts/utils/io_metadata.sh -m quality_control_raw_reads -s 0f24b004c_95a8cd511 -p PE
