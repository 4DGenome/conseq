time_start=$(date +"%s")
run_date=`date +"%Y-%m-%d-%H-%M"`
job_name=$pipeline_name-$pipeline_version

CONSEQ=/users/GR/mb/jquilez/projects/conseq

# python script to write/access metadata
io_metadata=$CONSEQ/scripts/utils/io_metadata.sh

# get species and assembly version from the metadata
species=homo_sapiens
if [[ $integrate_metadata == "yes" ]]; then
	sequencing_type=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a 'SEQUENCING_TYPE'`
	read_length=`$io_metadata -m get_from_metadata -s $sample_id -t input_metadata -a 'READ_LENGTH'`
	if [[ ${species,,} == 'homo_sapiens' ]]; then
		version=hg38_mmtv
	elif [[ ${species,,} == 'mus_musculus' ]]; then
		version=mm10
	fi
fi



#==================================================================================================
# PATHS
#==================================================================================================

# Directories

# pipeline scripts
SCRIPTS=$PIPELINE/scripts

# Primary output directory
if [[ $io_mode == "custom" ]]; then
	SAMPLE=$CUSTOM_OUT/$sample_id
else
	SAMPLE=$CONSEQ/data/$data_type/samples/$sample_id
fi

# Logs
LOGS=$SAMPLE/logs/$version

# Trim reads
TRIMMED=$SAMPLE/fastqs_processed/trimmomatic
SINGLE=$TRIMMED/single_end
PAIRED=$TRIMMED/paired_end
UNPAIRED=$TRIMMED/unpaired_reads
ADAPTERS=/software/mb/el6.3/Trimmomatic-0.33/adapters

# Mapping/alignment
GENOME_DIR=/users/GR/mb/jquilez/assemblies/$species/$version/star_genome_index/read_length_${read_length}bp
STAR=$SAMPLE/alignments/star/$version

# reads per transcript quantification
KALLISTO_QUANT=$SAMPLE/quantifications/kallisto/$version
FEATURECOUNTS_QUANT=$SAMPLE/quantifications/featurecounts/$version

# Reads per million (RPM) profiles
PROFILES=$SAMPLE/profiles/$version

# SHA cheksums
CHECKSUMS=$SAMPLE/checksums/$version/$run_date
checksums=$CHECKSUMS/files_checksums.sha


# Files

# input FASTQ
if [[ $io_mode == "custom" ]]; then
	ifq1_name=`grep $sample_id $CUSTOM_IN/sample_to_fastqs.txt |cut -f2`
	ifq2_name=`grep $sample_id $CUSTOM_IN/sample_to_fastqs.txt |cut -f3`
	ifq1=$CUSTOM_IN/$ifq1_name
	ifq2=$CUSTOM_IN/$ifq2_name
else
	ifq1=$CONSEQ/data/$data_type/raw/*/${sample_id}*read1.fastq.gz
	ifq2=$CONSEQ/data/$data_type/raw/*/${sample_id}*read2.fastq.gz
fi

# tools
java=`which java`
trimmomatic=`which trimmomatic`
star=`which star`
qualimap=`which qualimap`
samtools=`which samtools`
kallisto=`which kallisto`
featureCounts=`which featureCounts`
perl=`which perl`
bedGraphToBigWig=`which bedGraphToBigWig`

# indices and annotation
chrom_sizes=/users/GR/mb/jquilez/assemblies/$species/$version/ucsc/${version}_chr1-22XYMUn.chrom.sizes
if [[ $version == "hg19" || $version == "hg19_mmtv" ]]; then
	kallisto_index=/users/GR/mb/jquilez/assemblies/$species/hg19/kallisto_index/kallisto_${species}_hg19_ensGene.index
	transcripts_gtf=/users/GR/mb/jquilez/assemblies/$species/hg19/gencode/gencode.v19.annotation.gtf
elif [[ $version == "hg38" || $version == "hg38_mmtv" ]]; then
	kallisto_index=/users/GR/mb/jquilez/assemblies/$species/hg38/kallisto_index/kallisto_${species}_hg38_gencode_v24.index
	transcripts_gtf=/users/GR/mb/jquilez/assemblies/$species/hg38/gencode/gencode.v24.annotation.gtf
fi

# python script to write/access metadata
#io_metadata=/users/GR/mb/jquilez/utils/io_metadata.sh



# =================================================================================================
# CODE EXECUTION
# =================================================================================================

main() {

	echo 
	# store general parameters into the metadata
	if [[ $integrate_metadata == "yes" ]]; then
		$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a PIPELINE_RUN_MODE -v $pipeline_run_mode
		$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a QUEUE -v $queue
		$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a MEMORY -v $memory
		$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a MAX_TIME -v $max_time
		$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a SLOTS -v $slots
		$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a ASSEMBLY_VERSION -v $version
		$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a JOB_NAME -v $job_name		
		$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a PATH_JOB_FILE -v $path_job_file		
	fi

	if [[ $pipeline_run_mode == 'full' ]]; then
		trim_reads_trimmomatic
		align_star
		quality_alignments
		quantification_featurecounts
		quantification_kallisto
		clean_up
	elif [[ $pipeline_run_mode == 'trim_reads_trimmomatic' ]]; then trim_reads_trimmomatic
	elif [[ $pipeline_run_mode == 'align_star' ]]; then align_star
	elif [[ $pipeline_run_mode == 'quality_alignments' ]]; then quality_alignments
	elif [[ $pipeline_run_mode == 'quantification_featurecounts' ]]; then quantification_featurecounts
	elif [[ $pipeline_run_mode == 'quantification_kallisto' ]]; then quantification_kallisto
	elif [[ $pipeline_run_mode == 'clean_up' ]]; then clean_up
	fi
	echo

	# Final message
	message_info "pipeline" "completed successfully"
	message_time_pipeline 	

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
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a $field_name -v $length
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
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a $field_name -v $length
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
			step_log=$SAMPLE/logs/${sample_id}_${step}_single_end.log
			single1=$SINGLE/${sample_id}_read1.fastq.gz
			params="$ifq1 $single1"
			ODIR=$SINGLE
			if [[ $integrate_metadata == "yes" ]]; then
				$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a PATH_FASTQ_READ1 -v $ifq1
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
			step_log=$SAMPLE/logs/${sample_id}_${step}_paired_end.log
			paired1=$PAIRED/${sample_id}_read1.fastq.gz
			paired2=$PAIRED/${sample_id}_read2.fastq.gz
			unpaired1=$UNPAIRED/${sample_id}_read1.fastq.gz
			unpaired2=$UNPAIRED/${sample_id}_read2.fastq.gz
			params="$ifq1 $ifq2 $paired1 $unpaired1 $paired2 $unpaired2"
			ODIR=$PAIRED
			if [[ $integrate_metadata == "yes" ]]; then
				$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a PATH_FASTQ_READ1 -v $ifq1
				$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a PATH_FASTQ_READ2 -v $ifq2
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
	java -jar $trimmomatic $sequencing_type \
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
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a ADAPTERS_SEQS -v $seqs
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a SEED_MISMATCHES -v $seedMismatches
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a PALINDROME_CLIP_THRESHOLD -v $palindromeClipThreshold
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a SIMPLE_CLIP_THRESHOLD -v $simpleClipThreshold
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a MIN_ADAPTER_LENGTH -v $minAdapterLength
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a KEEP_BOTH_READS -v $keepBothReads
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a LEADING -v $leading
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a TRAILING -v $trailing
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a TARGET_LENGTH -v $targetLength
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a STRICTNESS -v $strictness
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a MIN_LENGTH -v $minLength
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_READS_TRIMMED -v $n_reads_trimmed
		message_info $step "trimmomatic parameters and numbe of trimmed reads added to metadata"
	fi

	# delete intermediate files
	message_info $step "trimmed reads are in $ODIR"
	if [[ $sequencing_type == "PE" ]]; then
		message_info $step "unpaired reads are deleted"
		rm -fr $UNPAIRED
	fi

	# data integrity
	mkdir -p $CHECKSUMS
	if [[ $sequencing_type == "SE" ]]; then
		shasum $ifq1 >> $checksums
	elif [[ $sequencing_type == "PE" ]]; then
		shasum $ifq1 >> $checksums
		shasum $ifq2 >> $checksums
	fi

	message_time_step $step $time0

}


# =================================================================================================
# Align reads with STAR
# =================================================================================================

align_star() {

	step="align_star"
	time0=$(date +"%s")

	# align paired-end reads with STAR
	star_version=`$star --version`
	message_info $step "align trimmed single-end reads with STAR (version = $star_version)"
	message_info $step "using ENCODE standard options for long RNA-seq pipeline"
	# Mapping options adjusted following ENCODE standard options for long RNA-seq found in the 
	# STAR manual https://github.com/alexdobin/STAR/tree/master/doc (note that although this links to the STAR 2.4 manual,
	# the options used here are also available for STAR 2.3, the version used here)
	# --outFilterType BySJout = reduces the number of ”spurious” junctions
	# --outFilterMultimapNmax 20 = max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
	# --alignSJoverhangMin 8 = minimum overhang for unannotated junctions
	# --alignSJDBoverhangMin 1 = minimum overhang for annotated junctions
	# --outFilterMismatchNmax 999 = maximum number of mismatches per pair, large number switches off this filter
	# --outFilterMismatchNoverLmax 0.04 = max number of mismatches per pair relative to read length: for 2x50b, 
	# max number of mis- matches is 0.04*100=4 for the pair (default is 0.3 so we are much more restrictive)
	# --alignIntronMin 20: ...
	# --alignIntronMax 1000000: maximum intron length
	# --outSAMtype BAM SortedByCoordinate = output sorted by coordinate
	# --outWigType: make read per million profile files. 4 files in total:
	# forward and reverse
	# unique alignments and multiple-alignments (provided the latter have less than outFilterMultimapNmax placements)  
	if [[ $sequencing_type == "SE" ]]; then
		step_log=$LOGS/${sample_id}_${step}_single_end.log
		single1=$SINGLE/${sample_id}_read1.fastq.gz
		ODIR1=$STAR/single_end
		ODIR2=$PROFILES/single_end
		params="$single1"
	elif [[ $sequencing_type == "PE" ]]; then
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
		paired1=$PAIRED/${sample_id}_read1.fastq.gz
		paired2=$PAIRED/${sample_id}_read2.fastq.gz
		ODIR1=$STAR/paired_end
		ODIR2=$PROFILES/paired_end
		params="$paired1 $paired2"
	fi
	mkdir -p $ODIR1
	mkdir -p $ODIR2
	TMP_DIR=$ODIR1/my_tmp
	$star \
		    --genomeDir $GENOME_DIR/ \
			--genomeLoad NoSharedMemory \
			--runThreadN $slots \
			--outFilterType "BySJout" \
			--outFilterMultimapNmax 20 \
			--alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 \
			--outFilterMismatchNmax 999 \
			--outFilterMismatchNoverLmax 0.04 \
			--alignIntronMin 20 \
			--alignIntronMax 1000000 \
			--alignMatesGapMax 1000000 \
			--readFilesIn $params \
			--outSAMtype BAM SortedByCoordinate \
			--outTmpDir $TMP_DIR/ \
			--outFileNamePrefix $ODIR1/$sample_id. \
			--outWigType bedGraph \
			--readFilesCommand zcat >$step_log 2>&1
	rm -fr $TMP_DIR
	message_info $step "alignments are in $ODIR1"

	# index BAM
	rm -f $ODIR1/$sample_id.Aligned.sortedByCoord.out.bam.bai
	$samtools index $ODIR1/$sample_id.Aligned.sortedByCoord.out.bam

	# convert bedGraph to bigWig (more suitable for UCSC browser upload)
	# unique alignments, strand1
	ibg=$ODIR1/$sample_id.Signal.Unique.str1.out.bg
	obw=$ODIR2/${sample_id}_unique_strand1_rpm.bw
	$bedGraphToBigWig $ibg $chrom_sizes $obw
	rm -f $ibg
	# unique alignments, strand2
	ibg=$ODIR1/$sample_id.Signal.Unique.str2.out.bg
	obw=$ODIR2/${sample_id}_unique_strand2_rpm.bw
	$bedGraphToBigWig $ibg $chrom_sizes $obw 
	rm -f $ibg
	# unique and multi alignments, strand1
	ibg=$ODIR1/$sample_id.Signal.UniqueMultiple.str1.out.bg
	obw=$ODIR2/${sample_id}_unique_multiple_strand1_rpm.bw
	$bedGraphToBigWig $ibg $chrom_sizes $obw
	rm -f $ibg
	# unique alignments, strand2
	ibg=$ODIR1/$sample_id.Signal.UniqueMultiple.str2.out.bg
	obw=$ODIR2/${sample_id}_unique_multiple_strand2_rpm.bw
	$bedGraphToBigWig $ibg $chrom_sizes $obw
	rm -f $ibg

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
	# parse
	star_log_final=$ODIR1/*Log.final.out
	# uniquely mapped reads
	n_reads_unique=`grep "Uniquely mapped reads number" $star_log_final |cut -f2 -d'|' | sed "s/\t//g"`
	p_reads_unique=`grep "Uniquely mapped reads %" $star_log_final |cut -f2 -d'|' | sed "s/%//g" | sed "s/\t//g"`
	# multi-mappings accepted (i.e. with less locations than specified with --outFilterMultimapNmax)
	n_reads_multi_mapping_accepted=`grep "Number of reads mapped to multiple loci" $star_log_final |cut -f2 -d'|'| sed "s/\t//g"`
	p_reads_multi_mapping_accepted=`grep "% of reads mapped to multiple loci" $star_log_final |cut -f2 -d'|' | sed "s/%//g"| sed "s/\t//g"`
	# multi-mappings excluded (i.e. with more locations than specified with --outFilterMultimapNmax)
	p_reads_multi_mapping_too_many=`grep "% of reads mapped to too many loci" $star_log_final |cut -f2 -d'|' | sed "s/%//g"| sed "s/\t//g"`
	# unmapped too short
	p_reads_unmapped_too_short=`grep "% of reads unmapped: too short" $star_log_final |cut -f2 -d'|' | sed "s/%//g"| sed "s/\t//g"`
	# number of splices detected
	n_splices=`grep "Number of splices: Total" $star_log_final |cut -f2 -d'|'| sed "s/\t//g"`
	message_info $step "reads unique (number) = $n_reads_unique"
	message_info $step "reads unique (percentage) = $p_reads_unique"
	message_info $step "accepted multi-mappings (number) = $n_reads_multi_mapping_accepted"
	message_info $step "accepted multi-mappings (percentage) = $p_reads_multi_mapping_accepted"
	message_info $step "excluded multi-mappings (percentage) = $p_reads_multi_mapping_too_many" 
	message_info $step "reads unmapped too short (percentage) = $p_reads_unmapped_too_short"
	message_info $step "splices (number) = $n_splices"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a GENOME_DIR -v $GENOME_DIR 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a OUT_FILTER_TYPE -v BySJout 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a OUT_FILTER_MULTIMAP_N_MAX -v 20 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a ALIGN_SJ_OVERHANG_MIN -v 8 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a ALIGN_SJDB_OVERHANG_MIN -v 1 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a OUT_FILTER_MISMATCH_N_MAX -v 999 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a OUT_FILTER_MISMATCH_N_OVER_L_MAX -v 0.04	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a ALIGN_INTRON_MIN -v 20 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a ALIGN_INTRON_MAX -v 1000000 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a ALIGN_MATES_GAP_MAX -v 1000000 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a VERSION_STAR -v $star_version 	
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_READS_UNIQUE -v $n_reads_unique 
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_READS_UNIQUE -v $p_reads_unique 
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_READS_MULTI_MAPPING_ACCEPTED -v $n_reads_multi_mapping_accepted 
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_READS_MULTI_MAPPING_ACCEPTED -v $p_reads_multi_mapping_accepted 
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_READS_MULTI_MAPPING_TOO_MANY -v $p_reads_multi_mapping_too_many 
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_READS_UNMAPPED_TOO_SHORT -v $p_reads_unmapped_too_short
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_SPLICES -v $n_splices
	fi

	# data integrity
	mkdir -p $CHECKSUMS
	shasum $ODIR1/$sample_id.*bam >> $checksums

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
		ibam=$STAR/single_end/$sample_id.Aligned.sortedByCoord.out.bam
		BAMQC=$STAR/single_end/qualimap_bamqc
		RNASEQ=$STAR/single_end/qualimap_rnaseq
	elif [[ $sequencing_type == "PE" ]]; then
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
		ibam=$STAR/paired_end/$sample_id.Aligned.sortedByCoord.out.bam
		BAMQC=$STAR/paired_end/qualimap_bamqc
		RNASEQ=$STAR/paired_end/qualimap_rnaseq
	fi

	# general QC of the BAM (bamqc)
	if [[ $strand_specific == 0 ]]; then p="non-strand-specific"
	elif [[ $strand_specific == 1 ]]; then p="forward-stranded"
	elif [[ $strand_specific == 2 ]]; then p="reverse-stranded"
	fi
	message_info $step "general QC of the BAM (using qualimap's bamqc)"
	$qualimap bamqc --java-mem-size=$memory -bam $ibam -outdir $BAMQC -c -ip -nt $slots -p $p >$step_log 2>&1

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
	genome_results=$BAMQC/genome_results.txt
	# globals
	n_mapped_paired_reads=`grep "number of mapped paired reads (first in pair)" $genome_results |cut -f2 -d"=" | sed "s/[ ,]//g"`
	n_overlapping_read_pairs=`grep "number of overlapping read pairs" $genome_results |cut -f2 -d"=" | sed "s/[ ,]//g"`
	p_overlapping_read_pairs=`echo "(100 * $n_overlapping_read_pairs) / $n_mapped_paired_reads" | bc -l`
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
	message_info $step "percentage of overlapping read pais = $p_overlapping_read_pairs"
	message_info $step "percentage duplication = $p_duplication"
	message_info $step "median insert size (bp) = $median_insert_size"
	message_info $step "mean mapping quality = $mean_mapping_quality"
	message_info $step "mean coverage (X) = $mean_coverage"
	message_info $step "mean coverage adjusted paired-end = $mean_coverage_paired_end"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_MAPPED_PAIRED_READS -v $n_mapped_paired_reads
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_OVERLAPPING_READ_PAIRS -v $p_overlapping_read_pairs
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_DUPLICATION -v $p_duplication
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a MEDIAN_INSERT_SIZE -v $median_insert_size
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a MEAN_MAPPING_QUALITY -v $mean_mapping_quality
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a MEAN_COVERAGE -v $mean_coverage
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a MEAN_COVERAGE_PAIRED_END -v $mean_coverage_paired_end
	fi

	# RNAseq-specific QC
	if [[ $strand_specific == 0 ]]; then p="non-strand-specific"
	elif [[ $strand_specific == 1 ]]; then p="strand-specific-forward"
	elif [[ $strand_specific == 2 ]]; then p="strand-specific-reverse"
	fi
	# sequencing type
	if [[ $sequencing_type == "PE" ]]; then
		pe="-pe"
	else
		pe=""
	fi
	echo
	message_info $step "RNAseq-specific QC of the BAM (using qualimap's rnaseq)"
	$qualimap rnaseq --java-mem-size=$memory -a proportional -bam $ibam -outdir $RNASEQ -gtf $transcripts_gtf -p $p $pe -s >$step_log 2>&1

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
	rnaseq_qc_results=$RNASEQ/rnaseq_qc_results.txt
	# reads alignment
	alignments_total=`grep "total alignments" $rnaseq_qc_results | cut -f2 -d'=' | sed "s/[ ,]//g"`
	alignments_secondary=`grep "secondary alignments" $rnaseq_qc_results | cut -f2 -d'=' | sed "s/[ ,]//g"`
	alignments_non_unique=`grep "non-unique alignments" $rnaseq_qc_results | cut -f2 -d'=' | sed "s/[ ,]//g"`
	alignments_to_genes=`grep "aligned to genes" $rnaseq_qc_results | cut -f2 -d'=' | sed "s/[ ,]//g"`
	alignments_no_feature_assigned=`grep "no feature assigned" $rnaseq_qc_results | cut -f2 -d'=' | sed "s/[ ,]//g"`
	not_aligned=`grep "not aligned" $rnaseq_qc_results | cut -f2 -d'=' | sed "s/[ ,]//g"`
	p_alignments_secondary=`echo "(100 * $alignments_secondary) / $alignments_total" | bc -l`
	p_alignments_non_unique=`echo "(100 * $alignments_non_unique) / $alignments_total" | bc -l`
	p_alignments_to_genes=`echo "(100 * $alignments_to_genes) / $alignments_total" | bc -l`
	p_alignments_no_feature_assigned=`echo "(100 * $alignments_no_feature_assigned) / $alignments_total" | bc -l`
	p_not_aligned=`echo "(100 * $not_aligned) / $alignments_total" | bc -l`
	# reads genomic origin
	p_alignments_exonic=`grep "exonic" $rnaseq_qc_results | cut -f2 -d'(' | sed "s/[,%)]//g"`
	p_alignments_intronic=`grep "intronic" $rnaseq_qc_results | cut -f2 -d'(' | sed "s/[,%)]//g"`
	p_alignments_intergenic=`grep "intergenic" $rnaseq_qc_results | cut -f2 -d'(' | sed "s/[,%)]//g"`
	p_alignments_overlapping_exon=`grep "overlapping exon" $rnaseq_qc_results | cut -f2 -d'(' | sed "s/[,%)]//g"`
	# print values
	message_info $step "total alignments = $alignments_total"
	message_info $step "percentage secondary alignments = $p_alignments_secondary"
	message_info $step "percentage non-unique alignments = $p_alignments_non_unique"
	message_info $step "percentage aligned to genes = $p_alignments_to_genes"
	message_info $step "percentage no feature assigned = $p_alignments_no_feature_assigned"
	message_info $step "percentage not aligned = $p_not_aligned"
	message_info $step "percentage exonic = $p_alignments_exonic"
	message_info $step "percentage intronic = $p_alignments_intronic"
	message_info $step "percentage intergenic = $p_alignments_intergenic"
	message_info $step "percentage overlapping exon = $p_alignments_overlapping_exon"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_TOTAL_ALIGNMENTS -v $alignments_total
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_ALIGNMENTS_SECONDARY -v $p_alignments_secondary
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_ALIGNMENTS_NON_UNIQUE -v $p_alignments_non_unique
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_ALIGNMENTS_TO_GENES -v $p_alignments_to_genes
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_ALIGNMENTS_NO_FEATURE_ASSIGNED -v $p_alignments_no_feature_assigned
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_NOT_ALIGNED -v $p_not_aligned
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_ALIGNMENTS_EXONIC -v $p_alignments_exonic
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_ALIGNMENTS_INTRONIC -v $p_alignments_intronic
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_ALIGNMENTS_INTERGENIC -v $p_alignments_intergenic
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_ALIGNMENTS_OVERLAPPING_EXON -v $p_alignments_overlapping_exon
	fi

	# remove log file for this step as this is very big (~40GB!)
	rm -f $step_log

	message_time_step $step $time0

}


# =================================================================================================
# Transcript quantification with featureCounts
# =================================================================================================

quantification_featurecounts() {

	step="quantification_featurecounts"
	time0=$(date +"%s")

	# quantify abundances of transcripts
	# -p = fragments will be counted instead of reads (the 2 paired-end reads are originated from 1 fragment)
	# -g gene_id = summarize reads counts per transcript
	# -s = strandness: 0=unstranded, 1=strand-specific, 2=reversely-stranded
	# ibam = alignments
	# transcripts_gtf = gene models
	message_info $step "quantifying read counts per gene using featureCounts"
	message_info $step "using gene models from $transcripts_gtf"
	message_info $step "kallisto (transcripts) and featurecounts (genes) quantifications are not directly comparable"
	message_info $step "even if the *_mmtv version of the assembly is used, gene/transcript models are the same as for *"
	if [[ $sequencing_type == "SE" ]]; then
		params=""
		IDIR=$STAR/single_end
		ODIR=$FEATURECOUNTS_QUANT/single_end
		message_info $step "sequencing type is $sequencing_type so reads are counted"
		step_log=$LOGS/${sample_id}_${step}_single_end.log
	elif [[ $sequencing_type == "PE" ]]; then
		params="-p"
		IDIR=$STAR/paired_end
		ODIR=$FEATURECOUNTS_QUANT/paired_end
		message_info $step "sequencing type is $sequencing_type so fragments are counted"
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
	fi
	ibam=$IDIR/$sample_id.Aligned.sortedByCoord.out.bam
	mkdir -p $ODIR
	ofile=$ODIR/${sample_id}_featurecounts.txt
	$featureCounts $params \
					-g gene_id \
					-s $strand_specific \
					-T $slots \
					-a $transcripts_gtf \
					-o $ofile \
					$ibam >$step_log 2>&1

	message_info $step "quantifications are in $ODIR"

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
	# parse
	log_final=$ODIR/*txt.summary
	fragments_total=`grep "Total fragments" $step_log | cut -f2 -d':' | sed "s/[ |]//g"`
	fragments_assigned=`grep "Assigned" $log_final | cut -f2`
	fragments_ambiguous=`grep "Ambiguity" $log_final | cut -f2`
	fragments_multimapping=`grep "MultiMapping" $log_final | cut -f2`
	fragments_no_features=`grep "NoFeatures" $log_final | cut -f2`
	p_fragments_assigned=`echo "(100 * $fragments_assigned) / $fragments_total" | bc -l`
	p_fragments_ambiguous=`echo "(100 * $fragments_ambiguous) / $fragments_total" | bc -l`
	p_fragments_multimapping=`echo "(100 * $fragments_multimapping) / $fragments_total" | bc -l`
	p_fragments_no_features=`echo "(100 * $fragments_no_features) / $fragments_total" | bc -l`
	message_info $step "total fragments = $fragments_total"
	message_info $step "percentage fragments assigned = $p_fragments_assigned"
	message_info $step "percentage fragments ambiguous = $p_fragments_ambiguous"
	message_info $step "percentage fragments multi-mapping = $p_fragments_multimapping"
	message_info $step "percentage fragments no features = $p_fragments_no_features"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_FRAGMENTS_TOTAL -v $fragments_total
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_FRAGMENTS_ASSIGNED -v $p_fragments_assigned
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_FRAGMENTS_AMBIGUOUS -v $p_fragments_ambiguous
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_FRAGMENTS_MULTIMAPPING -v $p_fragments_multimapping
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_FRAGMENTS_NO_FEATURES -v $p_fragments_no_features
	fi

	message_time_step $step $time0

}


# =================================================================================================
# Pseudo alignment and transcript quantification with kallisto
# =================================================================================================

quantification_kallisto() {

	step="quantification_kallisto"
	time0=$(date +"%s")

	# perform pseudoalignment and quantify abundances of transcripts
	# Using the option '--bias' produces an error
	# this has been reported and it should be fixed in the future:
	# https://groups.google.com/forum/#!searchin/kallisto-sleuth-users/bias|sort:date/kallisto-sleuth-users/d8nBIxrgESs/GSMKwhhfCgAJ
	kallisto_version=`$kallisto version`
	message_info $step "performing pseudoalignment and quantifying abundances of transcripts using $kallisto_version"
	message_info $step "using $kallisto_index as transcriptome reference"
	message_info $step "kallisto (transcripts) and featurecounts (genes) quantifications are not directly comparable"
	message_info $step "even if the *_mmtv version of the assembly is used, gene/transcript models are the same as for *"
	message_info $step "sequence based bias correction is only applied to single-end data, as it fails for paired-end"
	if [[ $sequencing_type == "SE" ]]; then
		ODIR=$KALLISTO_QUANT/single_end
		single1=$SINGLE/${sample_id}_read1.fastq.gz
		params="--single -l $fragment_length_avg -s $fragment_length_sd $single1 --bias"
		step_log=$LOGS/${sample_id}_${step}_single_end.log
		message_info $step "for single-end data, the user-provided fragment length average ($fragment_length_avg bp) is used"
		message_info $step "for single-end data, the user-provided fragment length standard deviation ($fragment_length_sd bp) is used"
	elif [[ $sequencing_type == "PE" ]]; then
		ODIR=$KALLISTO_QUANT/paired_end
		paired1=$PAIRED/${sample_id}_read1.fastq.gz
		paired2=$PAIRED/${sample_id}_read2.fastq.gz
		params="$paired1 $paired2 --bias"
		step_log=$LOGS/${sample_id}_${step}_paired_end.log
		message_info $step "for paired-end data, the fragment length average and standard deviation are inferred from the data"
	fi
	BOOTSTRAPS=$ODIR/bootstraps
	mkdir -p $BOOTSTRAPS
	$kallisto quant -i $kallisto_index -o $ODIR -t $slots -b $n_bootstraps $params >$step_log 2>&1
	message_info $step "converting form HDF5 to text"
	$kallisto h5dump -o $ODIR $ODIR/abundance.h5 >>$step_log 2>&1
	mv $ODIR/bs_abundance* $BOOTSTRAPS
	message_info $step "quantifications are in $ODIR"

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
	n_transcripts_target=`grep "number of targets:" $step_log | grep -v ',' |cut -f2 -d':' |sed "s/ //g"`
	n_reads_processed=`grep "pseudoaligned" $step_log |sed "s/ /;/g" |cut -f3 -d';' |sed "s/,//g"`
	n_reads_pseudoaligned=`grep "pseudoaligned" $step_log |sed "s/ /;/g" |cut -f5 -d';' |sed "s/,//g"`
	estimated_average_fragment_length=`grep "estimated average fragment length" $step_log | cut -f2 -d':' | sed "s/[ ,]//g"`
	p_reads_pseudoaligned=`echo "(100 * $n_reads_pseudoaligned) / $n_reads_processed" | bc -l`
	message_info $step "number transcripts quantified = $n_transcripts_target"
	message_info $step "number reads processed = $n_reads_processed"
	message_info $step "percentage reads pseudoaligned = $p_reads_pseudoaligned"
	message_info $step "estimated average fragment length (bp) = $estimated_average_fragment_length"

	# update metadata
	if [[ $integrate_metadata == "yes" ]]; then
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_TRANSCRIPTS_TARGET -v $n_transcripts_target
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_READS_PROCESSED -v $n_reads_processed
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a P_READS_PSEUDOALIGNED -v $p_reads_pseudoaligned
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a ESTIMATED_AVERAGE_FRAGMENT_LENGTH -v $estimated_average_fragment_length
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a KALLISTO_VERSION -v $kallisto_version
	 	$io_metadata -m add_to_metadata -t 'rnaseq' -s $sample_id -u $run_date -a N_BOOTSTRAPS -v $n_bootstraps
	fi

	message_time_step $step $time0

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
