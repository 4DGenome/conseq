# Managing the analysis of high-throughput sequencing data

Didactic dataset associated to the manuscript entitled "Managing the analysis of high-throughput sequencing data". The manuscript, cover letter, tables, figures (and more) can be found [here](https://drive.google.com/drive/folders/0B-MXr-KyKmm6VndUOXBqaGY5Vkk?usp=sharing).



<br>

## Table of Contents
- [Installation and usage](#installation-and-usage)
- [Didactic dataset](#didactic-dataset)
- [Sequencing index concordance](#sequencing-index-concordance)
- [Metadata collection](#metadata-collection)
- [Sample identifier (ID) scheme](#sample-identifier-(ID)-scheme)
- [Structured and hierarchical data organisation](#structured-and-hierarchical-data-organisation)
- [Scalability, parallelization and automatic configuration](#scalability-parallelization-and-automatic-configuration)
- [Dependencies](#dependencies)


<br>

## Installation and usage

Download the entire repository with:
```
git clone https://github.com/4DGenome/conseq.git
```

In many of the scripts in `scripts` paths are relative to the `$CONSEQ` Unix variable defined at the beginning of the script, which is set to `/users/GR/mb/jquilez/projects/conseq` (the absolute path to the repository directory in the machine where it was developed). As an example see:
```
head -n 12 scripts/utils/check_sequencing_index_concordance.sh
```

The scripts are written so that they can be executed from the directory where the repository is cloned by conveniently changing the `$CONSEQ` value, which can be achieved for all scripts with:
```
TARGET_DIR=my_home_directory
for s in scripts/*/*.sh; do
	IDIR='\/users\/GR\/mb\/jquilez\/projects\/conseq'
	sed -i "s/$IDIR/$TARGET_DIR/g" $s
done
```


<br>

## Didactic dataset

The Didactic dataset includes 7 high-throughput sequencing (HTS) samples (2 RNA-seq, 1 ChIP-seq and 4 Hi-C samples).

### HTS data

The `data` folder contains:
- FASTQ files with 1,000 and 2x1,000 reads for single and paired end samples, respectively, organised by HTS application and sequencing run date (`data/<data_type>/raw/<sequencing_run_date>`)
- data processed with core analysis pipelines (e.g. RNA-seq pipeline) (`data/<data_type>/samples`)

For information on how data are structured see [Structured and hierarchical data organisation](#structured-and-hierarchical-data-organisation).

### Metadata

The `metadata` directory contains the metadata for the 7 HTS samples in several formats:
- `metadata.db`: SQL database with several tables storing metadata collected for each HTS experiment as well as information resulting from running the analysis pipeline
- `metadata.tsv`: metadata as a tab-separated values file
- `metadata.xlsx`: metadata as a Microsoft Excel Open XML Format Spreadsheet file


<br>

## Sequencing index concordance

```bash
scripts/utils/check_sequencing_index_concordance.sh gv_066_01_01_chipseq
```
Executing the code above checks, for the sample passed as argument (`gv_066_01_01_chipseq in this case), whether the sequencing indexes found in the metadata and in the FASTQ agree. Note that for this sanity check to work:
1. a specific organisation of the data is required (see [Structured and hierarchical data organisation](#structured-and-hierarchical-data-organisation))
2. the sequencing index used in this sample needs to be part of the collected metadata (see [Metadata collection](#metadata-collection))
3. the sequencing index used in this sample needs to be part of the FASTQ header rows (which is not always provided); for instance:
```bash
zcat data/chipseq/raw/2012-11-29/gv_066_01_01_chipseq_read1.fastq.gz |head -n 4
```
Outputs:
```
@HWI-ST227:231:D1F2LACXX:1:1101:1228:2249 1:N:0:AACT
GGAGCTTTATTGAGTGTTAGAACAGCTCAGAGGAGATCCACAGTCA
+
FFFFHGHHHJJJJJIIIIJJJJJJJJJIJJJJJIIIJJJJJJJIJI
```
The first row shows that the sequencing index for this sample is AACT.



<br>

## Metadata collection

We provide scripts to download the collected metadata and dump them into a SQL database as well as to extract and update its content. In addition, the scripts allow printing dated freezes of the SQL database to monitor changes over time.


### Download input metadata 

```
scripts/utils/io_metadata.sh -m download_input_metadata
```
The `-m` command selects the mode. When `download_input_metadata` is passed, the metadata is downloaded from the corresponding URL and added to the `input_metadata` table of the metatadata SQL database.


### Extract metadata

```
scripts/utils/io_metadata.sh -m get_from_metadata -s gv_066_01_01_chipseq -t input_metadata -a SAMPLE_NAME
```
In addition to `-m get_from_metadata`, the following variables are required to extract a specific metadata value:
- sample (`-s`)
- table in the metadata database (`-t`)
- attribute (`-a`) that is printed out


### Update metadata

```
scripts/utils/io_metadata.sh -m add_to_metadata -s gv_066_01_01_chipseq -t input_metadata -a SAMPLE_NAME
```
Similarly, `-m add_to_metadata` updates a specific metadata value.

### Print freeze

```
scripts/utils/io_metadata.sh -m print_freeze
```
Prints the content of the SQL metadata database (one `*.csv` file per table) into a dated directory:
```
ls -l metadata/freezes/
```


<br>

## Sample identifier (ID) scheme

The SAMPLE_ID field in the SQL metadata database contains a unique sample ID.


### Manually assigned ID

In one of the two proposed sample ID schemes, the SAMPLE_ID value is manually assigned based on the metadata and contains 5 underscore-separated parts:
1. user (two-letter code corresponding to the first letter of the name and last name)
2. experiment ID (auto-incremental)
3. biological replicated ID (auto-incremental)
4. technical replicate ID (auto-incremental)
5. HTS application (e.g. rnaseq, chipseq)

For instance:
```
ifile=`ls metadata/freezes/*/input_metadata.csv`
grep 'APPLICATION\|fd' $ifile |column -t -s ',' |less -S
```
Selects the metadata for 2 biological replicates of a RNA-seq experiment performed by François Le Dily (fd).


### Automatically generated ID

In the other proposed sample ID scheme, the SAMPLE_ID value is automatically generated from the values in a set of specific fields. For instance:
```
scripts/utils/sample_id_generator.sh
```
Generates:
```
[info]	SAMPLE_ID for sample with Timestamp = 06/21/2016 17:21:57 is: 0f24b004c_f04ec1d87
[info]	SAMPLE_ID for sample with Timestamp = 09/21/2015 11:45:55 is: 0f24b004c_16d664eaa
[info]	SAMPLE_ID for sample with Timestamp = 11/09/2015 17:28:50 is: 66950b082_20b3c0316
[info]	SAMPLE_ID for sample with Timestamp = 11/09/2015 17:32:06 is: ad1a9f5b0_20b3c0316
```

From the metadata:
```
ifile=`ls metadata/freezes/*/input_metadata.csv`
grep 'APPLICATION\|4DGENOME' $ifile |column -t -s ',' |less -S
```
- the 2 first samples are 2 biological replicates and thus share the first half of the SAMPLE_ID
- the 2 last samples are samples from the same library and sequencing run


<br>

## Structured and hierarchical data organisation

We suggest a structured and hierarchical organisation that reflects the way in which sequencing data are generated and analyzed.


### Raw data

```
tree data/*/raw
```
As it can be seen from the output of the command above, FASTQ files (`*fastq.gz`) with the raw sequencing reads are named with the SAMPLE_ID and grouped by the run in which they were generated. Sequencing run directories contain not only the FASTQ files but also [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports informing about the quality of the raw reads.

### Processed data

```
tree data/*/samples
```

### Analysis results

```
tree projects
```
The command above shows the analysis performed for this paper, which is assigned to the `jquilez` project. The name of the analysis directory (`2017-04-07_analyses_manuscript`) starts with a timestamp plus a descriptive tag, so that additional analyses are naturally sorted chronologically. The analysis directory includes well-defined subdirectories where the output of the analysis is saved (`data`, `figures`, `scripts`, and `tables`). The analysis is thoroughly documented in the [Jupyter](http://jupyter.org/) notebook `2017-04-07_analyses_manuscript.ipynb`.


<br>

## Scalability, parallelization, automatic configuration and modularity

The diagram below illustrated incorporate scalability, parallelization, automatic configuration and modularity to our core analysis pipelines.

![Pipelines](https://github.com/CRG-Beato/utils/blob/master/pipelines.png)

Analysis pipelines can be simultaneously run on multiple samples with a single command (`*submit.sh *.config`). The configuration file (‘*.config’) contains the list of samples to be processed as well as the hard-coded parameters shared by all samples (e.g. number of processors or genome assembly version). The submission script (‘*submit.sh’) wraps the pipeline code (‘*seq.sh’) and generates, for each sample, a pipeline script with sample-specific variable values obtained from the SQL metadata database, and this will be submitted to the queuing system of the computing cluster where it will be executed (green). Selected metadata generated by the pipeline (e.g. running time, number of aligned reads) will be recorded into the database. For more flexibility, the pipeline modules to be executed are specified in the configuration file.

Specific details about the pipelines can be found at:
```
scripts/pipelines
```

<br>

## Dependencies
- wget
- Python and its packages:
	- os
	- sys
	- dataset
	- pandas
	- numpy
	- collections
	- glob
	- time
	- datetime
	- hashlib
- FastQC
- unzip



## To do
- [ ] update to zenodo:
	- `metadata/metadata.tsv`
	- `projects/jquilez/analysis/2017-04-07_analyses_manuscript/data/sra_stat_datasets.csv`
- [ ] make that metadata is effectively downloaded from Zenodo instead of importing the local file in `scripts/utils/io_metadata.sh`
- [ ] make that metadata is effectively downloaded from Zenodo instead of importing the local file in `scripts/utils/sample_id_generator.sh`
- [ ] add link to BSS17 slides
- [ ] add link to bioRxiv pre print


