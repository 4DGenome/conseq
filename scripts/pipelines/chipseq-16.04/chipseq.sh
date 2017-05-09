module load libpng

# additional run variables
time_start=$(date +"%s")
run_date=`date +"%Y-%m-%d-%H-%M"`
job_name=$pipeline_name-$pipeline_version

# python script to write/access metadata
io_metadata=/users/GR/mb/jquilez/utils/io_metadata.sh

# get species and assembly version from the metadata
species=homo_sapiens
if [[ $integrate_metadata == "yes" ]]; then
	sequencing_type=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a 'SEQUENCING_TYPE'`
	read_length=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a 'SEQUENCING_READ_LENGTH'`
	species=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a 'SPECIES'`	
	if [[ ${species,,} == 'homo_sapiens' ]]; then
		version=hg38_mmtv
	elif [[ ${species,,} == 'mus_musculus' ]]; then
		version=mm10
	fi
fi



#==================================================================================================
# PATHS
#==================================================================================================


# (1) Directories

# pipeline scripts
SCRIPTS=$PIPELINE/scripts

# Primary output directory
if [[ $io_mode == "custom" ]]; then
	SAMPLE=$CUSTOM_OUT/$sample_id
else
	SAMPLE=/users/GR/mb/jquilez/data/$data_type/samples/$sample_id
fi

# Logs
LOGS=$SAMPLE/logs/$version

# Trim reads
TRIMMED=$SAMPLE/fastqs_processed/trimmomatic
SINGLE=$TRIMMED/single_end
PAIRED=$TRIMMED/paired_end
UNPAIRED=$TRIMMED/unpaired_reads
ADAPTERS=/software/mb/el7.2/Trimmomatic-0.36/adapters

# alignment
BWA=$SAMPLE/alignments/bwa/$version

# tag directory
TAG_DIR=$SAMPLE/tag_directory/homer/$version

# reads per million (RPM) profiles
PROFILES=$SAMPLE/profiles/$version

# peaks
PEAKS=$SAMPLE/peaks/$peak_caller/$version

# SHA cheksums
CHECKSUMS=$SAMPLE/checksums/$version/$run_date
checksums=$CHECKSUMS/files_checksums.sha


# (2) Files

# input FASTQ
if [[ $io_mode == "custom" ]]; then
	ifq1_name=`grep -w $sample_id $CUSTOM_IN/sample_to_fastqs.txt |cut -f2`
	ifq2_name=`grep -w $sample_id $CUSTOM_IN/sample_to_fastqs.txt |cut -f3`
	ifq1=$CUSTOM_IN/$ifq1_name
	ifq2=$CUSTOM_IN/$ifq2_name
else
	ifq1=/users/GR/mb/jquilez/data/$data_type/raw/*/${sample_id}*read1.fastq.gz
	ifq2=/users/GR/mb/jquilez/data/$data_type/raw/*/${sample_id}*read2.fastq.gz
fi

# tools
trimmomatic=`which trimmomatic`
bwa=`which bwa`
qualimap=`which qualimap`
samtools=`which samtools`
makeTagDirectory=`which makeTagDirectory`
bedtools=`which bedtools`
perl=`which perl`
java=`which java`
bam2wig=`which bam2wig.pl`
bedGraphToBigWig=`which bedGraphToBigWig`
python=`which python`
bedgraph_to_bigwig=`which bedGraphToBigWig`

# genome fasta and chromosome sizes
if [[ ${species,,} == 'homo_sapiens' ]]; then
	genome_fasta=/users/GR/mb/jquilez/assemblies/${species,,}/$version/ucsc/${version}_chr1-22XYMUn.fa
	genome_chrom_sizes=/users/GR/mb/jquilez/assemblies/${species,,}/$version/ucsc/${version}_chr1-22XYMUn.chrom.sizes
elif [[ ${species,,} == 'mus_musculus' ]]; then
	genome_fasta=/users/GR/mb/jquilez/assemblies/${species,,}/$version/ucsc/${version}_chr1-19XYMUn.fa
	genome_chrom_sizes=/users/GR/mb/jquilez/assemblies/${species,,}/$version/ucsc/${version}_chr1-19XYMUn.chrom.sizes
fi



# =================================================================================================
# CODE EXECUTION
# =================================================================================================

main() {

	echo
	# store general parameters into the metadata
	if [[ $integrate_metadata == "yes" ]]; then
		$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PIPELINE_RUN_MODE -v $pipeline_run_mode
		$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a QUEUE -v $queue
		$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MEMORY -v $memory
		$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MAX_TIME -v $max_time
		$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a SLOTS -v $slots
		$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a ASSEMBLY_VERSION -v $version
		$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a JOB_NAME -v $job_name		
		$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PATH_JOB_FILE -v $path_job_file		
	fi

	if [[ $pipeline_run_mode == 'full' ]]; then
		trim_reads_trimmomatic
		align_bwa
		quality_alignments
		make_tag_directory
		make_profiles
		call_peaks
		clean_up
	elif [[ $pipeline_run_mode == 'full_no_call_peaks' ]]; then
		trim_reads_trimmomatic
		align_bwa
		quality_alignments
		make_tag_directory
		make_profiles
		clean_up
	elif [[ $pipeline_run_mode == 'full_no_make_profiles' ]]; then
		trim_reads_trimmomatic
		align_bwa
		quality_alignments
		make_tag_directory
		call_peaks
		clean_up
	elif [[ $pipeline_run_mode == 'full_from_alignments' ]]; then
		make_tag_directory
		make_profiles
		call_peaks
	elif [[ $pipeline_run_mode == 'trim_reads_trimmomatic' ]]; then trim_reads_trimmomatic
	elif [[ $pipeline_run_mode == 'align_bwa' ]]; then align_bwa
	elif [[ $pipeline_run_mode == 'quality_alignments' ]]; then quality_alignments
	elif [[ $pipeline_run_mode == 'make_tag_directory' ]]; then make_tag_directory
	elif [[ $pipeline_run_mode == 'make_profiles' ]]; then make_profiles
	elif [[ $pipeline_run_mode == 'call_peaks' ]]; then call_peaks
	elif [[ $pipeline_run_mode == 'clean_up' ]]; then clean_up
	fi
	echo

}



#==================================================================================================
# FUNCTIONS DEFINITIONS
#==================================================================================================


# =================================================================================================
# Pipeline progress functions
# =================================================================================================

# Outputs a message about the task being done
message_info() {
	step_name=$1
	message=$2
	echo -e "INFO \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
}

# Outputs a message about the error found and exits
message_error() {
	step_name=$1
	message=$2
	echo -e "ERROR \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
	exit	
}

# Outputs a warning message about the task being done
message_warn() {
	step_name=$1
	message=$2
	echo -e "WARN \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
}

# Outputs a message with the time in seconds spent in a given step
message_time_step() {
	step_name=$1
	field_name="TIME_${step_name^^}"
	time0=$2
	time1=$(date +"%s")
	length=$(($time1-$time0))
	echo -e "TIME \t`date +"%Y-%m-%d %T"` \t[$step_name] \tstep time for completion (seconds) = $length"
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a $field_name -v $length
	fi
	echo
}

# Outputs the total time in seconds for the pipeline to run
message_time_pipeline() {
	field_name="TIME_COMPLETE_PIPELINE"
	time0=$time_start
	time1=$(date +"%s")
	length=$(($time1-$time0))
	echo -e "TIME \t`date +"%Y-%m-%d %T"` \t[pipeline] \ttotal time for completion (seconds) = $length"
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a $field_name -v $length
	fi
	echo
}


# =================================================================================================
# Trim adapter and low-quality ends
# =================================================================================================

trim_reads_trimmomatic() {

	step="trim_reads_trimmomatic"
	time0=$(date +"%s")

	# Check that FASTQ file exists, make/define input/output directories and files
	if [[ $sequencing_type == "SE" ]]; then
		if [ -f $ifq1 ]; then
			mkdir -p $SAMPLE
			mkdir -p $SINGLE
			mkdir -p $LOGS
			mkdir -p $CHECKSUMS
			shasum $ifq1 >> $checksums
			step_log=$SAMPLE/logs/${sample_id}_${step}_single_end.log
			single1=$SINGLE/${sample_id}_read1.fastq.gz
			params="$ifq1 $single1"
			ODIR=$SINGLE
			if [[ $integrate_metadata == "yes" ]]; then
				$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PATH_FASTQ_READ1 -v $ifq1
				message_info $step "paths to read1 saved to metadata database"
			fi
		else
			message_error $step "$ifq1 not found. Exiting..."
		fi
	elif [[ $sequencing_type == "PE" ]]; then
		if [ -f $ifq1 ] && [ -f $ifq2 ]; then
			mkdir -p $SAMPLE
			mkdir -p $PAIRED
			mkdir -p $UNPAIRED
			mkdir -p $LOGS
			mkdir -p $CHECKSUMS
			shasum $ifq1 >> $checksums
			shasum $ifq2 >> $checksums
			step_log=$SAMPLE/logs/${sample_id}_${step}_paired_end.log
			paired1=$PAIRED/${sample_id}_read1.fastq.gz
			paired2=$PAIRED/${sample_id}_read2.fastq.gz
			unpaired1=$UNPAIRED/${sample_id}_read1.fastq.gz
			unpaired2=$UNPAIRED/${sample_id}_read2.fastq.gz
			params="$ifq1 $ifq2 $paired1 $unpaired1 $paired2 $unpaired2"
			ODIR=$PAIRED
			if [[ $integrate_metadata == "yes" ]]; then
				$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PATH_FASTQ_READ1 -v $ifq1
				$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PATH_FASTQ_READ2 -v $ifq2
				message_info $step "paths to read1 and read2 saved to metadata database"
			fi
		else
			message_error $step "$ifq1 and/or $ifq2 not found. Exiting..."
		fi
	fi

	# adapter trimming: the trimmomatic program directory contains a folder with the adapter sequences for
	# the Illumina sequencers in use. 'TruSeq3-PE.fa' is used, which contains the adapter sequences for the HiSeq
	message_info $step "sequencing type = $sequencing_type" 
	message_info $step "trimming adapter sequences for HiSeq, NextSeq or HiSeq"
	message_info $step "trimming low-quality reads ends using trimmomatic's recommended practices"
	seqs=$ADAPTERS/TruSeq3-$sequencing_type.fa
	targetLength=$read_length
	$java -jar $trimmomatic $sequencing_type \
					$params \
					ILLUMINACLIP:$seqs:$seedMismatches:$palindromeClipThreshold:$simpleClipThreshold:$minAdapterLength:$keepBothReads \
					LEADING:$leading \
					TRAILING:$trailing \
					MAXINFO:$targetLength:$strictness \
					MINLEN:$minLength >$step_log 2>&1

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
 	n_reads_trimmed=`grep Surviving $step_log | cut -f3 -d':' | cut -f1 -d'(' | sed "s/ //g"`
	message_info $step "reads after trimming = $n_reads_trimmed"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a ADAPTERS_SEQS -v $seqs
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a SEED_MISMATCHES -v $seedMismatches
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PALINDROME_CLIP_THRESHOLD -v $palindromeClipThreshold
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a SIMPLE_CLIP_THRESHOLD -v $simpleClipThreshold
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MIN_ADAPTER_LENGTH -v $minAdapterLength
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a KEEP_BOTH_READS -v $keepBothReads
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a LEADING -v $leading
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a TRAILING -v $trailing
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a TARGET_LENGTH -v $targetLength
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a STRICTNESS -v $strictness
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MIN_LENGTH -v $minLength
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a N_READS_TRIMMED -v $n_reads_trimmed
		message_info $step "trimmomatic parameters and numbe of trimmed reads added to metadata"
	fi

	# delete intermediate files
	message_info $step "trimmed reads are in $ODIR"
	if [[ $sequencing_type == "PE" ]]; then
		message_info $step "unpaired reads are deleted"
		rm -fr $UNPAIRED
	fi

	message_time_step $step $time0

}


# =================================================================================================
# aling reads with BWA
# =================================================================================================

align_bwa() {

	step="align_bwa"
	time0=$(date +"%s")

	# align reads with BWA
	message_info $step "align single-end reads with BWA"
	message_info $step "alignments are converted to BAM, sorted by genomic coordinates and multi-mappings filtered out"
	if [[ $sequencing_type == "SE" ]]; then
		step_log=$LOGS/${sample_id}_${step}_single_end.log
		single1=$SINGLE/${sample_id}_read1.fastq.gz
		ODIR=$BWA/single_end
		params="$single1"
	elif [[ $sequencing_type == "PE" ]]; then
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
		paired1=$PAIRED/${sample_id}_read1.fastq.gz
		paired2=$PAIRED/${sample_id}_read2.fastq.gz
		ODIR=$BWA/paired_end
		params="$paired1 $paired2"
	fi
	TMP_DIR=$BWA/my_tmp
	mkdir -p $TMP_DIR	
	mkdir -p $ODIR
	tbam=$ODIR/${sample_id}_sorted.bam
	obam=$ODIR/${sample_id}_sorted_filtered.bam
	read_group="@RG\tID:'$sample_id'\tLB:'$sample_id'\tPL:illumina\tPU:'$sample_id'\tSM:'$sample_id'"
	$bwa mem -t $slots -M $genome_fasta -R $read_group $params -v 0 |$samtools sort -o $tbam -O bam -T $TMP_DIR/$sample_id - >$step_log

	# clean reads:
	# unique mappings (implicit in -q 10 as BWA sets the mapping quality of multimappings to 0)
	# mapping quality >10 (-q 10)
	# exclude non primary alignments and supplementary alignments (-F 2304)
	# remove duplicates
	if [[ $sequencing_type == "SE" ]]; then
		$samtools view $tbam -bq 10 -F 2304 |$samtools rmdup -s - $obam >> $step_log
	elif [[ $sequencing_type == "PE" ]]; then
		$samtools view $tbam -bq 10 -F 2304 |$samtools rmdup - $obam >> $step_log
	fi
	$samtools index $obam

	# parse output
	n_alignments=`$samtools view $tbam | wc -l`
	n_alignments_filtered=`$samtools view $obam | wc -l`
	message_info $step "total number of alignments = $n_alignments"
	message_info $step "number of filtered alignments = $n_alignments_filtered"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a GENOME_FASTA -v $genome_fasta
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a N_ALIGNMENTS -v $n_alignments
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a N_ALIGNMENTS_FILTERED -v $n_alignments_filtered
		message_info $step "path to genome sequence FASTA and number of (uniquely) aligned reads added to metadata"
	fi

	# remove intermediate files
	rm -fr $TMP_DIR
	rm -f $tbam

	# data integrity
	mkdir -p $CHECKSUMS
	shasum $obam >> $checksums

	message_time_step $step $time0

}


# =================================================================================================
# Quality control of the mappings
# =================================================================================================

quality_alignments() {

	step="quality_alignments"
	time0=$(date +"%s")

	if [[ $sequencing_type == "SE" ]]; then
		step_log=$LOGS/${sample_id}_${step}_single_end.log
		ibam=$BWA/single_end/${sample_id}_sorted_filtered.bam
		BAMQC=$BWA/single_end/qualimap_bamqc
	elif [[ $sequencing_type == "PE" ]]; then
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
		ibam=$BWA/paired_end/${sample_id}_sorted_filtered.bam
		BAMQC=$BWA/paired_end/qualimap_bamqc
	fi

	# general QC of the BAM (bamqc)
	if [[ $strand_specific == 0 ]]; then p="non-strand-specific"
	elif [[ $strand_specific == 1 ]]; then p="forward-stranded"
	elif [[ $strand_specific == 2 ]]; then p="reverse-stranded"
	fi
	message_info $step "general QC of the BAM (using qualimap's bamqc)"
	$qualimap bamqc -bam $ibam -outdir $BAMQC -c -ip -nt $slots -p $p > $step_log 2>&1

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
	genome_results=$BAMQC/genome_results.txt
	# globals
	if [[ $sequencing_type == "PE" ]]; then
		n_mapped_paired_reads=`grep "number of mapped paired reads (first in pair)" $genome_results |cut -f2 -d"=" | sed "s/[ ,]//g"`
		n_overlapping_read_pairs=`grep "number of overlapping read pairs" $genome_results |cut -f2 -d"=" | sed "s/[ ,]//g"`
		p_overlapping_read_pairs=`echo "(100 * $n_overlapping_read_pairs) / $n_mapped_paired_reads" | bc -l`
		message_info $step "percentage of overlapping read pais = $p_overlapping_read_pairs"
	fi
	p_duplication=`grep "duplication rate" $genome_results | cut -f2 -d"=" | sed "s/[ %]//g"`
	# insert size
	median_insert_size=`grep "median insert size" $genome_results | cut -f2 -d"=" | sed "s/[ ,]//g"`
	# mapping quality
	mean_mapping_quality=`grep "mean mapping quality" $genome_results | cut -f2 -d"=" | sed "s/[ ,]//g"`
	# coverage (the paired-end coverage only makes sense for such kind of data)
	mean_coverage=`grep "mean coverageData" $genome_results | cut -f2 -d"=" | sed "s/[ X]//g"`
	if [[ $sequencing_type == "SE" ]]; then
		mean_coverage_paired_end='.'
	elif [[ $sequencing_type == "PE" ]]; then
		mean_coverage_paired_end=`grep "paired-end adapted mean coverage" $genome_results | cut -f2 -d"=" | sed "s/[ X]//g"`
	fi
	# print values
	message_info $step "percentage duplication = $p_duplication"
	message_info $step "median insert size (bp) = $median_insert_size"
	message_info $step "mean mapping quality = $mean_mapping_quality"
	message_info $step "mean coverage (X) = $mean_coverage"
	message_info $step "mean coverage adjusted paired-end = $mean_coverage_paired_end"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
		if [[ $sequencing_type == "PE" ]]; then
		 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a N_MAPPED_PAIRED_READS -v $n_mapped_paired_reads
		 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a P_OVERLAPPING_READ_PAIRS -v $p_overlapping_read_pairs
	 	fi
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a P_DUPLICATION -v $p_duplication
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MEDIAN_INSERT_SIZE -v $median_insert_size
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MEAN_MAPPING_QUALITY -v $mean_mapping_quality
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MEAN_COVERAGE -v $mean_coverage
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MEAN_COVERAGE_PAIRED_END -v $mean_coverage_paired_end
	fi

	# remove log file for this step as this is very big (~40GB!)
	rm -f $step_log

	message_time_step $step $time0

}


# =================================================================================================
# Make tag directory with HOMER
# =================================================================================================

make_tag_directory() {

	step="make_tag_directory"
	time0=$(date +"%s")

	# paths
	# this paths need to be created in case mode=full_from_mapping
	mkdir -p $LOGS
	mkdir -p $CHECKSUMS

	# to facilitate the analysis of ChIP-Seq (or any other type of short read re-sequencing data)
	# it is useful to first transform the sequence alignment into platform independent data structure representing the experiment
	# analogous to loading the data into a database.

	# Convert BAM to BED --required for making tag directory
	message_info $step "converting BAM to BED --required for making tag directory"
	if [[ $sequencing_type == "SE" ]]; then
		IDIR=$BWA/single_end
		ODIR=$TAG_DIR/single_end
		step_log=$LOGS/${sample_id}_${step}_single_end.log
	elif [[ $sequencing_type == "PE" ]]; then
		IDIR=$BWA/paired_end
		ODIR=$TAG_DIR/paired_end
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
	fi
	ibam=$IDIR/${sample_id}_sorted_filtered.bam
	obed=$IDIR/${sample_id}_sorted_filtered.bed
	$bedtools bamtobed -i $ibam | awk '{OFS="\t"; $4="."; print $0}' > $obed

	# Make tag directory with HOMER
	message_info $step "make tag directory with HOMER"
	mkdir -p $ODIR
	$makeTagDirectory $ODIR $obed -format bed -genome $genome_fasta -checkGC 2>$step_log

	# obed is only used by HOMER in the make_tag_directory step but this I file is no longer used
	# and of little use for the downstream analyses considering it is relatively big and can be
	# regenerated relatively quickly from the BAM if needed, so I will delete it
	rm -f $obed

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
	# parse
	estimated_genome_size=`grep "Estimated genome size" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	estimated_average_read_density_per_bp=`grep "Estimated average read density" $step_log | cut -f2 -d'=' |sed "s/ //g" |sed "s/perbp//g"`
	total_tags=`grep "Total Tags" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	total_positions=`grep "Total Positions" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	avg_tag_length=`grep "Average tag length" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	median_tags_per_position=`grep "Median tags per position" $step_log | cut -f2 -d'=' | cut -f1 -d'(' |sed "s/ //g"`
	avg_tags_per_position=`grep "Average tags per position" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	fragment_length_estimate=`grep "Fragment Length Estimate" $step_log | cut -f2 -d':' | sed "s/ //g"`
	peak_width_estimate=`grep "Peak Width Estimate" $step_log | cut -f2 -d':' | sed "s/ //g"`
	autocorrelation_same_strand_fold_enrichment=`grep "Same strand fold enrichment" $step_log | cut -f2 -d':' | sed "s/ //g"`
	autocorrelation_diff_strand_fold_enrichment=`grep "Diff strand fold enrichment" $step_log | cut -f2 -d':' | sed "s/ //g"`
	autocorrelation_same_to_diff_strand_fold_enrichment=`grep "Same / Diff fold enrichment" $step_log | cut -f2 -d':' | sed "s/ //g"`
	avg_fragment_gc=`grep "Avg Fragment GC" $step_log | cut -f2 -d'=' | sed "s/ //g" |sed "s/%//g"`
	avg_expected_gc=`grep "Avg Expected GC" $step_log | cut -f2 -d'=' | sed "s/ //g" |sed "s/%//g"`
	# print
	message_info $step "estimated genome size = $estimated_genome_size"
	message_info $step "estimated average read density per bp = $estimated_average_read_density_per_bp"
	message_info $step "total tags = $total_tags"
	message_info $step "total_positions = $total_positions"
	message_info $step "avg. tag length = $avg_tag_length"
	message_info $step "median tags per position = $median_tags_per_position"
	message_info $step "avg. tags per position = $avg_tags_per_position"
	message_info $step "fragment length estimate = $fragment_length_estimate"
	message_info $step "peak width estimate = $peak_width_estimate"
	message_info $step "autocorrelation: same strand fold enrichment = $autocorrelation_same_strand_fold_enrichment"
	message_info $step "autocorrelation: diff strand fold enrichment = $autocorrelation_diff_strand_fold_enrichment"
	message_info $step "autocorrelation: same-to-diff strand fold enrichment = $autocorrelation_same_to_diff_strand_fold_enrichment"
	message_info $step "avg. fragment GC% = $avg_fragment_gc"
	message_info $step "avg. expected GC% = $avg_expected_gc"
	
	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a ESTIMATED_GENOME_SIZE -v $estimated_genome_size
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a ESTIMETED_AVERAGE_READ_DENSITY_PER_BP -v $estimated_average_read_density_per_bp
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a TOTAL_TAGS -v $total_tags
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a TOTAL_POSITIONS -v $total_positions
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a AVG_TAG_LENGTH -v $avg_tag_length
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a MEDIAN_TAGS_PER_POSITION -v $median_tags_per_position
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a AVG_TAGS_PER_POSITION -v $avg_tags_per_position
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a FRAGMENT_LENGTH_ESTIMATE -v $fragment_length_estimate
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PEAK_WIDTH_ESTIMATE -v $peak_width_estimate
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a AUTOCORRELATION_SAME_STRAND_FOLD_ENRICHMENT -v $autocorrelation_same_strand_fold_enrichment
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a AUTOCORRELATION_DIFF_STRAND_FOLD_ENRICHMENT -v $autocorrelation_diff_strand_fold_enrichment
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a AUTOCORRELATION_SAME_TO_DIFF_STRAND_FOLD_ENRICHMENT -v $autocorrelation_same_to_diff_strand_fold_enrichment
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a AVG_FRAGMENT_GC -v $avg_fragment_gc
	 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a AVG_EXPECTED_GC -v $avg_expected_gc
		message_info $step "ChIP-seq metrics calculated with Homer added to metadata"
	fi

	message_time_step $step $time0

}



# =================================================================================================
# Make read per million (RPM) profiles
# =================================================================================================

make_profiles() {

	step="make_profiles"
	time0=$(date +"%s")

	# Generate RPM fragment profile
	message_info $step "generate reads per million profile (RPM) fragment profile"
	if [[ $sequencing_type == "SE" ]]; then
		# paths
		IDIR=$BWA/single_end
		ODIR=$PROFILES/single_end
		step_log=$LOGS/${sample_id}_${step}_single_end.log
		make_tag_directory_log=$LOGS/${sample_id}_make_tag_directory_single_end.log
		tag_info=$TAG_DIR/single_end/tagInfo.txt	
		ibam=$IDIR/${sample_id}_sorted_filtered.bam
		mkdir -p $ODIR
		# get the fragment length estimate
		fragment_length_estimate=`grep "Fragment Length Estimate" $make_tag_directory_log | cut -f2 -d':' | sed "s/ //g"`
		fragment_length_estimate_corrected=`cat $tag_info | grep fragmentLengthEstimate |cut -f 2 -d"=" | sed 's/[^0-9]*//g'`
		# get the number of million mapped reads
		n_mapped=`$samtools idxstats $ibam |cut -f3 |paste -sd+ |bc`
		scale_factor=`$python -c "print(1000000. / $n_mapped)"`
		# calculate the coverage (in bedGraph format)
		tbg=$ODIR/tmp.bgedGraph
		$bedtools genomecov -bg -ibam $ibam -g $genome_chrom_sizes -scale $scale_factor > $tbg
		# convert bedGraph to bigwig
		orpm=$ODIR/$sample_id.rpm.bw
		$bedgraph_to_bigwig $tbg $genome_chrom_sizes $orpm > $step_log
		# remove intermediate reads
		rm -f $tbg

	elif [[ $sequencing_type == "PE" ]]; then
		# paths
		IDIR=$BWA/paired_end
		ODIR=$PROFILES/paired_end
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
		make_tag_directory_log=$LOGS/${sample_id}_make_tag_directory_paired_end.log
		tag_info=$TAG_DIR/paired_end/tagInfo.txt	
		ibam=$IDIR/${sample_id}_sorted_filtered.bam
		mkdir -p $ODIR
		# get the fragment length estimate
		fragment_length_estimate=`grep "Fragment Length Estimate" $make_tag_directory_log | cut -f2 -d':' | sed "s/ //g"`
		fragment_length_estimate_corrected=`cat $tag_info | grep fragmentLengthEstimate |cut -f 2 -d"=" | sed 's/[^0-9]*//g'`
		# get valid pairs
		tbam=$ODIR/tmp.bam
 		$samtools view -bf 0x2 $ibam > $tbam
 		$samtools index $tbam
		# get the number of million mapped reads
		n_mapped=`$samtools idxstats $tbam |cut -f3 |paste -sd+ |bc`
		scale_factor=`$python -c "print(1000000. / $n_mapped)"`
		# calculate the coverage (in bedGraph format)
		tbg=$ODIR/tmp.bgedGraph
		$bedtools genomecov -bg -ibam $tbam -g $genome_chrom_sizes -scale $scale_factor > $tbg
		# convert bedGraph to bigwig
		orpm=$ODIR/$sample_id.rpm.bw
		$bedgraph_to_bigwig $tbg $genome_chrom_sizes $orpm > $step_log
 		rm $tbam $tbam.bai $tbg

	fi

	# data integrity
	mkdir -p $CHECKSUMS
	shasum $orpm >> $checksums

	message_time_step $step $time0

}



# =================================================================================================
# Call peaks
# =================================================================================================

call_peaks() {

	step="call_peaks"
	time0=$(date +"%s")

	if [[ $sequencing_type == "SE" ]]; then
		IDIR=$BWA/single_end
		step_log=$LOGS/${sample_id}_${step}_single_end.log
		tag_info=$TAG_DIR/single_end/tagInfo.txt	
	elif [[ $sequencing_type == "PE" ]]; then
		IDIR=$BWA/paired_end
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
		tag_info=$TAG_DIR/paired_end/tagInfo.txt	
	fi
	ibam=$IDIR/${sample_id}_sorted_filtered.bam

	# Get fragment length estimate as calculated in the 'make_tag_directory' step
	message_info $step "Get fragment length estimate (l) as calculated in the 'make_tag_directory' step"
	fragment_length_estimate_corrected=`cat $tag_info | grep fragmentLengthEstimate |cut -f 2 -d"=" | sed 's/[^0-9]*//g'`
	message_info $step "Fragment length (l) is $fragment_length_estimate_corrected bp (note this is not used if peak caller is zerone)"

	# Peak calling with MACS2
	if [[ $peak_caller == "macs2" ]]; then
		
		macs2=`which $peak_caller`
		message_info $step "Peak calling with MACS2"
		if [[ ${species,,} == "homo_sapiens" ]]; then
			genome_size="hs"
		elif [[ ${species,,} == "mus_musculus" ]]; then
			genome_size="mm"
		fi
		message_info $step "genome size for $species will be used"
		message_info $step "q-value cutoff = $macs2_qvalue (default is 0.01)"
		message_info $step "--nomodel (MACS2 will not try to model peak length but use l = $fragment_length_estimate_corrected instead"
		message_info $step "--call-summits = MACS reanalyzes the shape of signal profile to deconvolve subpeaks within each peak called"

		# Peak calling with input DNA as control
		if [[ $use_control == "yes" ]]; then
			if [ ! -f $control_bam ]; then
				message_error $step "$control_bam not found. Exiting..."
			fi
			if [[ $sequencing_type == "SE" ]]; then
			 	ODIR=$PEAKS/with_control/single_end
			 	step_log=$LOGS/${sample_id}_${step}_${peak_caller}_with_control_single_end.log
			elif [[ $sequencing_type == "PE" ]]; then
			 	ODIR=$PEAKS/with_control/paired_end
			 	step_log=$LOGS/${sample_id}_${step}_${peak_caller}_with_control_paired_end.log
			fi		
		 	mkdir -p $ODIR
		 	message_info $step "peak calling with input DNA ($control_bam) as control"
		 	$macs2 callpeak -t $ibam \
		 					-c $control_bam \
		 					-q $macs2_qvalue \
		 					--nomodel \
		 					--extsize $fragment_length_estimate_corrected \
		 					--outdir $ODIR \
		 					--call-summits \
		 					--gsize $genome_size \
		 					--name $sample_id 2>$step_log

			# parse step log to extract generated metadata
	 		n_peaks=`cat $ODIR/*.narrowPeak | wc -l`
			message_info $step "peaks = $n_peaks"

	 		# update metadata
		 	if [[ $integrate_metadata == "yes" ]]; then
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PEAK_CALLER -v $peak_caller
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a QVALUE_CUTOFF -v $macs2_qvalue
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PATH_CONTROL_BAM -v $control_bam
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a FRAGMENT_LENGTH_ESTIMATE_CORRECTED -v $fragment_length_estimate_corrected
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a N_PEAKS -v $n_peaks
				message_info $step "peak calling with input DNA as control added to metadata"
			fi

			# data integrity
			mkdir -p $CHECKSUMS
			shasum $ODIR/*.narrowPeak >> $checksums

		# Peak calling with the sample alone
		elif [[ $use_control == "no" ]]; then
			if [[ $sequencing_type == "SE" ]]; then
			 	ODIR=$PEAKS/sample_alone/single_end
			 	step_log=$LOGS/${sample_id}_${step}_${peak_caller}_sample_alone_single_end.log
			elif [[ $sequencing_type == "PE" ]]; then
			 	ODIR=$PEAKS/sample_alone/paired_end
			 	step_log=$LOGS/${sample_id}_${step}_${peak_caller}_sample_alone_paired_end.log
			fi		
			mkdir -p $ODIR
			message_info $step "peak calling with the sample alone (i.e. no input)"
			$macs2 callpeak -t $ibam \
							-q $macs2_qvalue \
							--nomodel \
							--extsize $fragment_length_estimate_corrected \
							--outdir $ODIR \
		 					--call-summits \
		 					--gsize $genome_size \
							--name $sample_id 2>$step_log
	
			# parse step log to extract generated metadata
	 		n_peaks=`cat $ODIR/*.narrowPeak | wc -l`
			message_info $step "peaks = $n_peaks"

	 		# update metadata
		 	if [[ $integrate_metadata == "yes" ]]; then
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PEAK_CALLER -v $peak_caller
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a QVALUE_CUTOFF -v $macs2_qvalue
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a FRAGMENT_LENGTH_ESTIMATE_CORRECTED -v $fragment_length_estimate_corrected
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a N_PEAKS -v $n_peaks
				message_info $step "peak calling without input DNA as control added to metadata"
			fi

			# data integrity
			mkdir -p $CHECKSUMS
			shasum $ODIR/*.narrowPeak >> $checksums

		fi

	# Peak calling with Zerone
	elif [[ $peak_caller == "zerone" ]]; then

		zerone=`which $peak_caller`
		message_info $step "peak calling with zerone"

		if [[ $use_control == "yes" ]]; then
			if [ ! -f $control_bam ]; then
				message_error $step "$control_bam not found. Exiting..."
			fi
			if [[ $sequencing_type == "SE" ]]; then
			 	ODIR=$PEAKS/with_control/single_end
			 	step_log=$LOGS/${sample_id}_${step}_${peak_caller}_with_control_single_end.log
			elif [[ $sequencing_type == "PE" ]]; then
			 	ODIR=$PEAKS/with_control/paired_end
			 	step_log=$LOGS/${sample_id}_${step}_${peak_caller}_with_control_paired_end.log
			fi		
			mkdir -p $ODIR
			message_info $step "peak calling with input DNA ($control_bam) as control"
			# Zerone parameters
			# -l = zerone produces an alternative output in which 
			# (i) only enriched windows (i.e. peaks) are shown
			# (ii) contiguous windows are merged
			otab1=$ODIR/${sample_id}_zerone.txt
			otab2=$ODIR/${sample_id}_zerone_enriched_merged.txt
			message_info $step "making discretization"
			$zerone -c -0 $control_bam -1 $ibam > $otab1
			message_info $step "making discretization only printing enriched regions"
			$zerone -l -0 $control_bam -1 $ibam > $otab2

			# parse step log to extract generated metadata
	 		n_peaks=`grep -v "#" $otab2 | wc -l`
			message_info $step "peaks = $n_peaks"

	 		# update metadata
		 	if [[ $integrate_metadata == "yes" ]]; then
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PEAK_CALLER -v $peak_caller
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a PATH_CONTROL_BAM -v $control_bam
			 	$io_metadata -m add_to_metadata -t 'chipseq' -s $sample_id -u $run_date -a N_PEAKS -v $n_peaks
				message_info $step "peak calling with input DNA as control added to metadata"
			fi

			# data integrity
			mkdir -p $CHECKSUMS
			shasum $otab1 >> $checksums
			shasum $otab2 >> $checksums

		# A control is mandatory for zerone, thus exiting otherwise
		else
			message_error $step "a control is mandatory for zerone. Exiting..."
		fi
	fi

}


# ========================================================================================
# Deletes intermediate files
# ========================================================================================

clean_up() {

	step="clean_up"
	time0=$(date +"%s")

	message_info $step "deleting the following intermediate files/directories:"
	message_info $step "$SAMPLE/fastqs_processed/trimmomatic/*/*"
	rm -f $SAMPLE/fastqs_processed/trimmomatic/*/*
	message_time_step $step $time0

}


main


