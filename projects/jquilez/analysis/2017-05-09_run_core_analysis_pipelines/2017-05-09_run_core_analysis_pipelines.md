# 2017-05-09_run_core_analysis_pipelines
----------------------------------------------------------------------------------------------------

**objective: run the RNA-seq and ChIP-seq pipelines**

**paths are relative to /users/GR/mb/jquilez/projects/conseq**



# [2017-05-09] Run RNA-seq pipeline
----------------------------------------------------------------------------------------------------

Configuration file used:
```
; This configuration file follows the INI file format (https://en.wikipedia.org/wiki/INI_file)

[data_type]
data_type			= rnaseq

[samples]
samples				=fd_005_01_01_rnaseq fd_005_02_01_rnaseq 				; e.g.: `samples=s1 s2 s3`, use 'test1' for testing purposes

[pipeline]
pipeline_name		= rnaseq
pipeline_version	= 16.06
pipeline_run_mode	= full

[IO mode]
io_mode				= standard									; standard = /users/GR/mb/jquilez, custom = in and out dir specified
CUSTOM_IN			= scripts/pipelines/rnaseq-16.06/test 		; only used if pipeline_io_mode=custom
CUSTOM_OUT			= scripts/pipelines/rnaseq-16.06/test		; only used if pipeline_io_mode=custom
sample_to_fastqs	= sample_to_fastqs.txt				; file with paths, relative to CUSTOM_IN, to read1 (and read2) FASTS, only used if pipeline_io_mode=custom

[cluster options]
submit_to_cluster	= no					; the following are only applied if submit_to_cluster=yes
queue				= long-sl7			; for real data = long-sl65, for test = short-sl65
memory				= 60G					; for real data = 60G, for test = 40G
max_time			= 24:00:00 				; for real data = 24:00:00, for test = 1:00:00
slots				= 8 					; for real data = 8, for test = 1
email				= javier.quilez@crg.eu	; email to which start/end/error emails are sent

[metadata]
integrate_metadata	= yes					; yes: metadata is stored into database

[trimmomatic]
; for recommended values see http://www.broadinstitute.org/cancer/software/genepattern/modules/docs/Trimmomatic/
; and those used in the supplementary data of the Trimmomatic paper (Bolger et al. 2014)
sequencing_type		= 					; PE=paired-end, SE=single-end
seedMismatches			= 2
palindromeClipThreshold	= 30
simpleClipThreshold		= 12
leading					= 3
trailing				= 3
minAdapterLength		= 1
keepBothReads			= true
minQual					= 3
strictness				= 0.999
minLength				= 36

[star]
species				= 
version				= 
read_length			= 

[kallisto]
n_bootstraps			= 100
fragment_length_avg		= 150				; for single-end data this is required; for paired-end data it is inferred from the data
fragment_length_sd		= 30				; for single-end data this is required; for paired-end data it is inferred from the data

[featureCounts]
strand_specific			= 2 				; 0=unstranded, 1=stranded, 2=reversely stranded
```

Execute pipeline:
```
scripts/pipelines/rnaseq-16.06/rnaseq_submit.sh scripts/pipelines/rnaseq-16.06/rnaseq.config
```



# [2017-05-09] Run ChIP-seq pipeline
----------------------------------------------------------------------------------------------------

Configuration file used:
```
; This configuration file follows the INI file format (https://en.wikipedia.org/wiki/INI_file)

[data_type]
data_type			= chipseq

[samples]
samples				=gv_066_01_01_chipseq			; e.g.: `samples=s1 s2 s3`, use 'test' for testing purposes

[pipeline]
pipeline_name		= chipseq
pipeline_version	= 16.04
pipeline_run_mode	= full

[IO mode]
io_mode				= standard									; standard = /users/GR/mb/jquilez, custom = in and out dir specified
CUSTOM_IN			= scripts/pipelines/chipseq-16.04/test 		; only used if pipeline_io_mode=custom
CUSTOM_OUT			= scripts/pipelines/chipseq-16.04/test 		; only used if pipeline_io_mode=custom
sample_to_fastqs	= sample_to_fastqs.txt				; file with paths, relative to CUSTOM_IN, to read1 (and read2) FASTS, only used if pipeline_io_mode=custom

[cluster options]
submit_to_cluster	= no					; the following are only applied if submit_to_cluster=yes
queue				= long-sl7				; for real data = long-sl65, for test = short-sl65
memory				= 80G					; for real data = 60G, for test = 40G
max_time			= 48:00:00 				; for real data = 48:00:00, for test = 1:00:00
slots				= 10 					; for real data = 10, for test = 1
email				= javier.quilez@crg.eu	; email to which start/end/error emails are sent

[metadata]
integrate_metadata	= yes					; yes: metadata is stored into database

[trimmomatic]
; for recommended values see http://www.broadinstitute.org/cancer/software/genepattern/modules/docs/Trimmomatic/
; and those used in the supplementary data of the Trimmomatic paper (Bolger et al. 2014)
sequencing_type		= 					; PE=paired-end, SE=single-end
seedMismatches			= 2
palindromeClipThreshold	= 30
simpleClipThreshold		= 12
leading					= 3
trailing				= 3
minAdapterLength		= 1
keepBothReads			= true
minQual					= 3
strictness				= 0.999
minLength				= 36
read_length				= 				; required if integrate_metadata=no, otherwise, ignored and retrieved from the metadata

[bwa]
species				= 			; required if integrate_metadata=no, otherwise, ignored and retrieved from the metadata
version				= 				; required if integrate_metadata=no, otherwise, ignored: hg38_mmtv (homo_sapiens) and mm10 (mus_musculus)

[qualimap]
strand_specific			= 1 				; 0=unstranded, 1=stranded, 2=reversely stranded

[peak calling]
peak_caller			= macs2
use_control			= no
control_bam			= 			; Input DNA to be used as control in the peak calling (e.g. $HOME/projects/er_pr/data/alignments/input_mcf7_merged_sorted_dups_removed.bam)

[macs2]
macs2_qvalue		= 0.05					; MACS2 defaults is 0.01 and 0.05 for broad marks
```


Execute pipeline:
```
scripts/pipelines/chipseq-16.04/chipseq_submit.sh scripts/pipelines/chipseq-16.04/chipseq.config
```

