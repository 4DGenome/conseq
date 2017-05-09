#!/bin/bash
#$ -N quality_control_fd_005_02_01_rnaseq
#$ -q short-sl7
#$ -l virtual_free=15G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/projects/conseq/projects/jquilez/analysis/2017-05-09_run_core_analysis_pipelines/job_out/quality_control_fd_005_02_01_rnaseq_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/conseq/projects/jquilez/analysis/2017-05-09_run_core_analysis_pipelines/job_out/quality_control_fd_005_02_01_rnaseq_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/conseq/data/rnaseq/raw/2014-06-06/fd_005_02_01_rnaseq_read1.fastq.gz -o /users/GR/mb/jquilez/projects/conseq/data/rnaseq/raw/2014-06-06/fastqc; rm -f /users/GR/mb/jquilez/projects/conseq/data/rnaseq/raw/2014-06-06/fastqc/fd_005_02_01_rnaseq*read1_fastqc.zip
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/conseq/data/rnaseq/raw/2014-06-06/fd_005_02_01_rnaseq_read2.fastq.gz -o /users/GR/mb/jquilez/projects/conseq/data/rnaseq/raw/2014-06-06/fastqc; rm -f /users/GR/mb/jquilez/projects/conseq/data/rnaseq/raw/2014-06-06/fastqc/fd_005_02_01_rnaseq*read2_fastqc.zip
/users/GR/mb/jquilez/projects/conseq/scripts/utils/io_metadata.sh -m quality_control_raw_reads -s fd_005_02_01_rnaseq -p PE
