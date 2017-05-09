# README


**Pipeline to call binding site peaks using ChIP-seq data**

<br>

## Modules

The pipeline is broken down into modules:

1. `trim_reads_trimmomatic`: trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. `align_bwa`: align reads with [BWA](http://bio-bwa.sourceforge.net/bwa.shtml)
3. `quality_alignments`: quality control of the mappings using [qualimap](http://qualimap.bioinfo.cipf.es/)
4. `make_tag_directory`: generate quality reads of the ChIP-seq fragments with [HOMER](http://homer.salk.edu/homer/)
5. `make_profiles`: generate reads per million (RPM) profiles
6. `call_peaks`: call binding site peaks with (a) [MACS2](https://github.com/taoliu/MACS) (either using control or not) or (b) [Zerone](https://github.com/gui11aume/zerone) (control is required)
7. `clean_up`: delete relatively large intermediate files which can be re-generated

The modules can be executed altogether or individually (see [Configuration file](#configuration-file)). The diagram below shows the order in which modules are sequentially executed (numbers), when the full pipeline is run, and the dependencies between modules in case they want to be run individually (e.g. all modules require that `trim_reads_trimmomatic` has been executed):

![chipseq-16.04](https://github.com/4DGenome/conseq/blob/master/docs/figures_github_repo/chipseq-16.04/chipseq-16.04.001.png)


<br>

## Scripts

- `chipseq.sh`: script with the code of the pipeline
- `chipseq_submit`: wrapper script that both:
	- retrieves configuration variables and parameter values from the `chipseq.config` file
	- (if applies) submits jobs (one per sample) to execute the pipeline in a Univa Grid Engine HPC cluster 
- `chipseq.config`: configuration file with the list of samples and the hard-coded parameter values (see [Configuration file](#configuration-file))


<br>

## Pipeline execution

```
chipseq_submit.sh chipseq.config
```


<br>

## Configuration file

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

- `pipeline_run_mode`:
	- by setting `trim_reads_trimmomatic`, `align_bwa`, `quality_alignments`, `make_tag_directory`, `make_profiles`, `call_peaks` or `clean_up` the corresponding module (see above) is executed (**note modules are sequential so each cannot be run unless the preceding one is completed**) 
	- `full`: all the steps above in the sequential order
	- `full_no_call_peaks`: all the modules except `call_peaks` (this is useful when a sample will be used as control)
	- `full_no_make_profiles`: all the modules except `make_profiles`
	- `full_from_alignments`: all modules after `quality_alignmnents`
- `data_type`: chipseq, atacseq (if `io_mode=standard` this will be used to search for the input/output directory)
- `samples`: space-separated list of samples
- `io_mode`:
	- `standard`: pre-defined in/out directories within the `$CONSEQ` path
	- `custom`:	in/out directories defined by the user in `CUSTOM_IN` and `CUSTOM_OUT`, and the input FASTQs defined in `sample_to_fastqs.txt`
- `[cluster options]`: submission options for an Univa Grid Engine HPC cluster (e.g. see [CRG cluster](http://www.linux.crg.es/index.php/Main_Page))
- `[trimmomatic]`: see [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- `species` and `version`: as of 2016-03-01, genome sequence FASTA files are available for `homo_sapiens`, `version` hg19 and hg38 (either with or without the MMTV construct)
- `peak_caller`: macs2, zerone
- `use_control`: yes/no to using a BAM as control or input when calling peaks
- `control_bam`: BAM file to be used as control or input; mandatory if `use_control=yes`


<br>

## `io_mode` and `integrate_metadata`

When `io_mode = standard`, FASTQ input file(s) and output directory are pre-defined. This behaviour can be changed with `io_mode = custom`, where the `chipseq.config` file provides the path to the input FASTQs and a file listing them as well as the path to the output directory. As this is a non-standard usage of the pipeline, `integrate_metadata` is internally set to **no** so values for all the variables in the `chipseq.config` are required. To try the custom mode one can use the data in the `test` directory.


<br>

## Dependencies and edits

- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- edit `$ADAPTERS` so that it points to the `adapters` subdirectory of the downloaded version of Trimmomatic
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml)
- [Java](https://www.java.com/en/)
- [QualiMap](http://qualimap.bioinfo.cipf.es/)
- [SAMtools](http://samtools.sourceforge.net/)
- [HOMER's](http://homer.ucsd.edu/homer/) `makeTagDirectory` tool
- [BEDtools](http://bedtools.readthedocs.io/en/latest/)
- [Perl](https://www.perl.org/)
- [kentUtils](https://github.com/ENCODE-DCC/kentUtils) (`bedGraphToBigWig`)
- [Python](https://www.python.org/)
- edit `genome_fasta` so that it points to the genome assembly reference sequence
- edit `chrom_sizes` so that it points to the file with the size of the chromosomes
