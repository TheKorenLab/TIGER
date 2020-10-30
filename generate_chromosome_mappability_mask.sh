#! /bin/bash

# Argument parser
while [ "$1" != "" ]; do
    case $1 in
        --ref )    shift
                   REF_GENOME=$1
                   ;;
       --k )       shift
                   KMER_LENGTH=$1
                   ;;
      --chr )      shift
                   CHR=$1
                   ;;
    esac
    shift
done

# Generate list of all possible reads
./simulate_chromosome_reads.py --ref $REF_GENOME'.fa' --k $KMER_LENGTH --chr $CHR

# Sort reads
java -jar $PICARD SortSam INPUT=$CHR'.ubam' OUTPUT=$CHR'_sorted.ubam' SORT_ORDER=queryname 2> /dev/null

# Mark adapters
java -jar $PICARD MarkIlluminaAdapters INPUT=$CHR'_sorted.ubam' OUTPUT=$CHR'_sorted_adapters.ubam' \
    METRICS=/dev/null 2> /dev/null

# Align synthetic reads to reference genome
java -jar $PICARD SamToFastq INPUT=$CHR'_sorted_adapters.ubam' FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 \
    INTERLEAVE=true NON_PF=true 2> /dev/null| \
    bwa mem -M -t 7 -p $REF_GENOME /dev/stdin 2> /dev/null | \
    java -jar $PICARD MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=$CHR'_sorted.ubam' \
        OUTPUT=$CHR'_sorted_piped.ubam' REFERENCE_SEQUENCE=$REF_GENOME'.fa' CREATE_INDEX=true ADD_MATE_CIGAR=true \
        CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS 2> /dev/null

# Mark duplicates
java -jar $PICARD MarkDuplicates INPUT=$CHR'_sorted_piped.ubam' OUTPUT=$CHR'.bam' METRICS_FILE=/dev/null \
    REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT 2> /dev/null

# Index BAM
java -jar $PICARD BuildBamIndex INPUT=$CHR'.bam' 2> /dev/null

# Output list of mappable positions
if [ ${CHR:0:3} = 'chr' ]; then
    OUTFILE='mappability_mask_'$CHR'.txt'
else
    OUTFILE='mappability_mask_chr'$CHR'.txt'
fi

samtools view -F4 -F16 -F1024 -q1  $CHR'.bam' $CHR | cut -f4 > $OUTFILE

# Clean up
rm $CHR* header.sam
