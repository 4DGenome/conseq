#!/bin/bash
#$ -N quality_control_gv_066_01_01_chipseq
#$ -q short-sl7
#$ -l virtual_free=15G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/projects/conseq/projects/jquilez/analysis/2017-05-09_run_core_analysis_pipelines/job_out/quality_control_gv_066_01_01_chipseq_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/projects/conseq/projects/jquilez/analysis/2017-05-09_run_core_analysis_pipelines/job_out/quality_control_gv_066_01_01_chipseq_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1
/software/mb/bin/fastqc --extract /users/GR/mb/jquilez/projects/conseq/data/chipseq/raw/2012-11-29/gv_066_01_01_chipseq_read1.fastq.gz -o /users/GR/mb/jquilez/projects/conseq/data/chipseq/raw/2012-11-29/fastqc; rm -f /users/GR/mb/jquilez/projects/conseq/data/chipseq/raw/2012-11-29/fastqc/gv_066_01_01_chipseq*read1_fastqc.zip
/users/GR/mb/jquilez/projects/conseq/scripts/utils/io_metadata.sh -m quality_control_raw_reads -s gv_066_01_01_chipseq -p SE
