#!/bin/bash


ml Trimmomatic
ml STAR
ml RSEM

echo "Please enter the absolute path to the input directory with the raw data files in fastq format:"
read INPUT_DIR
echo "Please enter the absolute path to the output directory:"
read OUTPUT_DIR
echo "Please enter the absolute path to the genome directory with the STAR indexes:"
read GENOME_DIR
echo "Please enter the absolute path to the RSEM reference directory: \
(Make sure to include the common prefix in the path \
e.g. if your RSEM reference was built with human(hg38) annotation then \
the path would be <path_to_rsem_ref>/hg38)"
read RSEM_DIR

mkdir -p $OUTPUT_DIR/fastqc
fastqc -t 5 -o $OUTPUT_DIR/fastqc/ $INPUT_DIR/*.gz

echo "Please enter the number of threads:"
read THREADS

mkdir -p $OUTPUT_DIR/trimmed_files

# set parameters for lines 7 and 8

A=$INPUT_DIR/hcc1395_normal_rep1_r1.fastq.gz  # Replace with the actual file names .fastq.gz or .fq.gz
B=$INPUT_DIR/hcc1395_normal_rep1_r2.fastq.gz 
C=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r1_paired.fq.gz # STAR uses paired files
D=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r1_unpaired.fq.gz 
E=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r2_paired.fq.gz # STAR uses paired files
F=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r2_unpaired.fq.gz 
F=ILLUMINACLIP:$INPUT_DIR/TruSeq3-PE.fa:2:30:10  # TruSeq3-PE.fa is the adapter file in form .fa or .fasta
G=LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $THREADS $A $B $C $D $E $F $G

mkdir -p $OUTPUT_DIR/trimmed_files/fastqc
fastqc $OUTPUT_DIR/trimmed_files/fastqc/hcc1395_normal_rep1_r1_paired.fq.gz
# do all files in $OUTPUT_DIR/trimmed_files


mkdir -p $OUTPUT_DIR/STAR_alignment

A=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r1_pa 
B=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r2_pa
C=--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN 10
D=--outFileNamePrefix $OUTPUT_DIR/STAR_alignment/hcc1395_normal_rep1 #change the prefix to match the sample name

STAR --genomeDir $GENOME_DIR --readFilesCommand zcat --readFilesIn $A $B $C $D

mkdir -p $OUTPUT_DIR/RSEM_counts

rsem-calculate-expression --bam --no-bam-output -p 5 --paired-end --forward-prob 0 \
$OUTPUT_DIR/STAR_alignment/hcc1395_normal_rep1Aligned.toTranscriptome.out.bam \
$RSEM_DIR  \ 
$OUTPUT_DIR/RSEM_counts/hcc1395_normal_rep1 >& $OUTPUT_DIR/RSEM_counts/hcc1395_normal_rep1.log # change the file name to match the sample