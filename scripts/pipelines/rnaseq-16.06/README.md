# README

**Pipeline to quantify gene/transcript abundance using RNA-seq data**


<br>

## Modules


The pipeline is broken down into modules:

1. `trim_reads_trimmomatic`: trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. `align_star`: align reads to genome with [STAR](https://github.com/alexdobin/STAR)
3. `quality_alignments`: quality control of the mappings using [qualimap](http://qualimap.bioinfo.cipf.es/)
4. `quantification_featurecounts`: quantification of reads counts per gene, using alignments from `align_star`, with [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
5. `quantification_kallisto`: pseudo-alignment of reads + quantification of read counts per transcript with [Kallisto](http://pachterlab.github.io/kallisto/)
6. `clean_up`: delete relatively large intermediate files which can be re-generated

The modules can be executed altogether or individually (see [Configuration file](##configuration-file)). The diagram below shows the order in which modules are sequentially executed (numbers), when the full pipeline is run, and the dependencies between modules in case they want to be run individually (e.g. all modules require that `trim_reads_trimmomatic` has been executed):

![rnaseq-16.04](https://github.com/4DGenome/conseq/blob/master/docs/figures_github_repo/rnaseq-16.04/rnaseq-16.04.001.png)


<br>

## Scripts

- `rnaseq.sh`: script with the code of the pipeline
- `rnaseq_submit`: wrapper script that both:
	- retrieves configuration variables and parameter values from the `rnaseq.config` file
	- (if applies) submits jobs (one per sample) to execute the pipeline in an Univa Grid Engine HPC cluster 
- `rnaseq.config`: configuration file with the list of samples and the hard-coded parameter values (see [Configuration file](##configuration-file))


<br>

## Pipeline execution

```
rnaseq_submit.sh rnaseq.config
```


<br>

## Configuration file


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

- `data_type`: `rnaseq`
- `samples`: espace-separated list of samples
- `pipeline_run_mode`: any of the 5 steps described in the **Modules** section or `full` to run all of them sequentially
- `io_mode`:
	- `standard`: pre-defined in/out directories within the `$CONSEQ` path
	- `custom`:	in/out directories defined by the user in `CUSTOM_IN` and `CUSTOM_OUT`, and the input FASTQs defined with `read1_fname` and, if paired-end, `read2_fname` (see below)
- `[cluster options]`: submission options for an Univa Grid Engine HPC cluster (e.g. see [CRG cluster](http://www.linux.crg.es/index.php/Main_Page))
- `[trimmomatic]`: see [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- `species` and `version`: as of 2016-03-01, STAR genome files are only generated for `homo_sapiens`, `version` hg19 and hg38, and `read_length` 50, 75 and 100 bp
- `[kallisto]`: see [Kallisto](http://pachterlab.github.io/kallisto/)
- `[featureCounts]`: see [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)


<br>

## `io_mode` and `integrate_metadata`

When `io_mode = standard`, FASTQ input file(s) and output directory are pre-defined. This behaviour can be changed with `io_mode = custom`, where the `rnaseq.config` file provides the path to the input FASTQs and a file listing them as well as the path to the output directory. As this is a non-standard usage of the pipeline, `integrate_metadata` is internally set to **no** so values for all the variables in the `rnaseq.config` are required. To try the custom mode one can use the data in the `test` directory.


<br>

## Dependencies and edits

- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- edit `$ADAPTERS` so that it points to the `adapters` subdirectory of the downloaded version of Trimmomatic
- [STAR](https://github.com/alexdobin/STAR)
- edit `$GENOME_DIR` so that it points to the STAR-indexed genome files
- [Java](https://www.java.com/en/)
- [QualiMap](http://qualimap.bioinfo.cipf.es/)
- [SAMtools](http://samtools.sourceforge.net/)
- [Kallisto](http://pachterlab.github.io/kallisto/)
- edit `kallisto_index` so that it points to the Kallisto-indexed transcriptome
- [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
- [Perl](https://www.perl.org/)
- [kentUtils](https://github.com/ENCODE-DCC/kentUtils) (bedGraphToBigWig)
- edit `transcripts_gtf` so that it points to the genome assembly annotation
- edit `chrom_sizes` so that it points to the file with the size of the chromosomes

*For the sake of comparison, use the transcriptome sequences used by Kallisto and the genome assembly annotation from the same source (e.g. [GENCODE](http://www.gencodegenes.org/)).*