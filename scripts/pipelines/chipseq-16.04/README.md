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

The modules can be executed altogether or individually (see [Configuration file](##configuration-file)). The diagram below shows the order in which modules are sequentially executed (numbers), when the full pipeline is run, and the dependencies between modules in case they want to be run individually (e.g. all modules require that `trim_reads_trimmomatic` has been executed):

![rnaseq-16.04](https://github.com/4DGenome/conseq/blob/master/docs/figures_github_repo/chipseq-16.04/chipseq-16.04.001.png)

## Scripts

- `chipseq.sh`: most of the code
- `chipseq_submit`: wrapper script that both:
	- retrieves configuration variables and parameter values from the `chipseq.config` file
	- (if applies) submits jobs (one per sample) to execute the pipeline in the CRG cluster
- `chipseq.config`: configuration file (see below)


## Execute pipeline

```
/pipeline_location/chipseq_submit.sh <*.config>
```

**Users other than me have no writting permissions for the `chipseq.config` file, so they need to provide their own file**


## Configuration file

- `pipeline_run_mode`:
	- by setting `trim_reads_trimmomatic`, `align_bwa`, `quality_alignments`, `make_tag_directory`, `make_profiles`, `call_peaks` or `clean_up` the corresponding module (see above) is executed (**note modules are sequential so each cannot be run unless the preceding one is completed**) 
	- `full`: all the steps above in the sequential order
	- `full_no_call_peaks`: all the modules except `call_peaks` (this is useful when a sample will be used as control)
	- `full_from_alignments`: all modules after `quality_alignmnents`
- `data_type`: chipseq, atacseq (if `io_mode=standard` this will be used to search for the input/output directory)
- `samples`: space-separated list of samples
- `pipeline_run_mode`: any of the 5 steps described in the **Modules** section of `full` to run all of them sequentially
- `io_mode`:
	- `standard`: in/out directories defined in my home directory `/users/GR/mb/jquilez`
	- `custom`:	in/out directories defined by the user in `CUSTOM_IN` and `CUSTOM_OUT`, and the input FASTQs defined in `sample_to_fastqs.txt`
- `[cluster options]`: see [CRG cluster](http://www.linux.crg.es/index.php/Main_Page)
- **if `integrate_metadata=yes` (see `*.config`), do not execute any two steps of the pipeline simulateneously as it will block the database**
- `[trimmomatic]`: see [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- `species` and `version`: as of 2016-03-01, genome sequence FASTA files are available for `homo_sapiens`, `version` hg19 and hg38 (either with or without the MMTV construct)
- `peak_caller`: macs2, zerone
- `use_control`: yes/no to using a BAM as control or input when calling peaks
- `control_bam`: BAM file to be used as control or input; mandatory if `use_control=yes`


## Test dataset

2 samples downloaded from the SRA data repository:
```
IDIR=data/atacseq/raw/2016-03-01
ODIR=pipelines/chipseq-16.03/test
# 1860 reads
cp -v $IDIR/SRR1779699_read1.fastq.gz $ODIR/test1_read1.fastq.gz
cp -v $IDIR/SRR1779699_read2.fastq.gz $ODIR/test1_read2.fastq.gz
# 1962 reads
cp -v $IDIR/SRR1779697_read1.fastq.gz $ODIR/test2_read1.fastq.gz
cp -v $IDIR/SRR1779697_read2.fastq.gz $ODIR/test2_read2.fastq.gz
```