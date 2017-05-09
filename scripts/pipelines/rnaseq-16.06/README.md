# README

**Pipeline to quantify gene/transcript abundance using RNA-seq data**


<br>

## Modules

The pipeline is broken down into modules that are run sequentially and that can be run individually too.

1. `trim_reads_trimmomatic`: trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. `align_star`: align reads to genome with [STAR](https://github.com/alexdobin/STAR)
3. `quality_alignments`: quality control of the mappings using [qualimap](http://qualimap.bioinfo.cipf.es/)
4. `quantification_featurecounts`: quantification of reads counts per gene, using alignments from `align_star`, with [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
5. `quantification_kallisto`: pseudo-alignment of reads + quantification of read counts per transcript with [Kallisto](http://pachterlab.github.io/kallisto/)
6. `clean_up`: delete relatively large intermediate files which can be re-generated

The diagram below shows the order in which modules are sequentially executed (numbers) and the dependencies between modules (e.g. all modules require that `trim_reads_trimmomatic` has been executed):

![rnaseq-16.04](https://github.com/4DGenome/conseq/blob/master/docs/figures_github_repo/rnaseq-16.04/rnaseq-16.04.001.png)


<br>

## Scripts

- `rnaseq.sh`: modules code
- `rnaseq_submit`: wrapper script that both:
	- retrieves configuration variables and parameter values from the `rnaseq.config` file
	- (if applies) submits jobs (one per sample) to execute the pipeline in an Univa Grid Engine HPC cluster 
- `rnaseq.config`: configuration file (see below)


<br>

## Pipeline execution

```
rnaseq_submit.sh rnaseq.config
```


<br>

## Configuration file

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