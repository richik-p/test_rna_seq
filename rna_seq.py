
import os
import subprocess
import argparse
import gzip
import tarfile

# Rest of the argparse setup

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description='RNA-Seq analysis pipeline using STAR for alignment.')
parser.add_argument('--data_dir', required=True, help='Directory containing your FASTQ files')
parser.add_argument('--output_dir', required=True, help='Output directory for results')
parser.add_argument('--ref_genome_dir', required=True, help='Directory of the STAR reference genome indices')
parser.add_argument('--rsem', required=True, help='Path to the RSEM reference directory')
parser.add_argument('--ann', required=True, help='Path to the annotation GTF file')
parser.add_argument('--unzip', action='store_true', help='Unzip .gz files before processing')
# Parse the command-line arguments
args = parser.parse_args()

# Function to unzip .gz files
def unzip_file(file_path, output_path):
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as gz_file:
            with open(output_path, 'w') as out_file:
                for line in gz_file:
                    out_file.write(line)
    elif file_path.endswith('.tar'):
        with tarfile.open(file_path, 'r') as tar_file:
            tar_file.extractall(path=output_path, filter='data')

# Assign arguments to variables
data_directory = args.data_dir
output_directory = args.output_dir
reference_genome_dir = args.ref_genome_dir
rsem_dir = args.rsem
annotation_file_path = args.ann
unzip = args.unzip

# Assuming args.unzip will be True if the --unzip flag is used
if args.unzip:
    # Modify the directory listing to include .gz files and unzip them
    fastq_files = [f for f in os.listdir(data_directory) if f.endswith('.fastq.gz') or f.endswith('.fq.gz')]
    for gz_file in fastq_files:
        # Define the full path to the compressed file
        gz_file_path = os.path.join(data_directory, gz_file)
        # Define the output file path (same name without .gz)
        output_file_path = os.path.join(data_directory, gz_file[:-3])
        # Unzip the file
        unzip_file(gz_file_path, output_file_path)
        # Replace the .gz file in the list with the uncompressed file
        fastq_files[fastq_files.index(gz_file)] = gz_file[:-3]
else:
    # Original listing for non-compressed files
    fastq_files = [f for f in os.listdir(data_directory) if f.endswith('.fastq') or f.endswith('.fq')]


# Define paths to your tools
star_path = "STAR"
stringtie_path = "stringtie"
samtools_path = "samtools"
fastqc = "fastqc"

# Define the output directory for your results
output_directory = os.path.join(os.getcwd(), "output_folder")
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Define the OUTPUT_DIR/fastqc directory for storing FastQC results if not already present
fastqc_output_dir = os.path.join(output_directory, "fastqc")
if not os.path.exists(fastqc_output_dir):
    os.makedirs(fastqc_output_dir)

# fastqc -t 5 -o $OUTPUT_DIR/fastqc/ $INPUT_DIR/*.gz in pyhton
subprocess.run([fastqc, "-t", "5", "-o", fastqc_output_dir, os.path.join(data_directory, "*.gz")])

# mkdir -p $OUTPUT_DIR/trimmed_files in python if not already present
trimmed_files_dir = os.path.join(output_directory, "trimmed_files")
if not os.path.exists(trimmed_files_dir):
    os.makedirs(trimmed_files_dir)

# Find all pairs of FASTQ files ending with r1.fastq.gz/r1.fq.gz and the names of files before r1/r2 and get n/2 names of files
fastq_files = [f for f in os.listdir(data_directory) if f.endswith('.fastq.gz') or f.endswith('.fq.gz')]
true_names_and_file = [(f.split('r1')[0], f) for f in fastq_files if f.endswith('r1.fastq.gz') or f.endswith('r1.fq.gz')]

#mkdir -p $OUTPUT_DIR/trimmed_files in python if not already present
trimmed_files_dir = os.path.join(output_directory, "trimmed_files")
if not os.path.exists(trimmed_files_dir):
    os.makedirs(trimmed_files_dir)


# find pairs of r1 and r2 files from fastq_files using true_names
for n in true_names_and_file:
    r1_file = n[1]
    r2_file = n[0] + 'r2' + r1_file.split('r1')[1]
    name = n[0].rstrip('._')

    subprocess.run(["java", "-jar", "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar", "PE", "-threads", "8", 
                    os.path.join(data_directory, r1_file), os.path.join(data_directory, r2_file), 
                    os.path.join(trimmed_files_dir, name + "_r1_paired.fq.gz"), 
                    os.path.join(trimmed_files_dir, name + "_r1_unpaired.fq.gz"), 
                    os.path.join(trimmed_files_dir, name + "_r2_paired.fq.gz"), 
                    os.path.join(trimmed_files_dir, name + "_r2_unpaired.fq.gz"), 
                    "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"])
    
    # mkdir -p $OUTPUT_DIR/trimmed_files/fastqc in python if not already present
    trimmed_fastqc_output_dir = os.path.join(trimmed_files_dir, "fastqc")
    if not os.path.exists(trimmed_fastqc_output_dir):
        os.makedirs(trimmed_fastqc_output_dir)

    # fastqc $OUTPUT_DIR/trimmed_files/fastqc/hcc1395_normal_rep1_r1_paired.fq.gz in python
    subprocess.run([fastqc, os.path.join(trimmed_files_dir, name + "_r1_paired.fq.gz")])

    # mkdir -p $OUTPUT_DIR/STAR_alignment in python if not already present
    alignment_output_dir = os.path.join(output_directory, name + "_STAR")
    if not os.path.exists(alignment_output_dir):
        os.makedirs(alignment_output_dir)

    # A=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r1_paired 
    # B=$OUTPUT_DIR/trimmed_files/hcc1395_normal_rep1_r2_paired
    # C=--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM --runThreadN 10
    # D=--outFileNamePrefix $OUTPUT_DIR/STAR_alignment/hcc1395_normal_rep1 #change the prefix to match the sample name

    # STAR --genomeDir $GENOME_DIR --readFilesCommand zcat --readFilesIn $A $B $C $D

    subprocess.run([star_path, "--genomeDir", reference_genome_dir, "--readFilesCommand", "zcat", 
                    "--readFilesIn", os.path.join(trimmed_files_dir, name + "_r1_paired.fq.gz"), 
                    os.path.join(trimmed_files_dir, name + "_r2_paired.fq.gz"), "--outSAMtype BAM SortedByCoordinate", 
                    "--outSAMunmapped Within", "--twopassMode Basic", "--outFilterMultimapNmax 1", "--quantMode TranscriptomeSAM",
                    "--runThreadN", "10", "--outFileNamePrefix", os.path.join(alignment_output_dir, name)])
    
    # mkdir -p $OUTPUT_DIR/RSEM_counts in python if not already present
    rsem_output_dir = os.path.join(output_directory, "RSEM_counts")
    if not os.path.exists(rsem_output_dir):
        os.makedirs(rsem_output_dir)

    # rsem-calculate-expression --bam --no-bam-output -p 5 --paired-end --forward-prob 0 \
    # $OUTPUT_DIR/STAR_alignment/hcc1395_normal_rep1Aligned.toTranscriptome.out.bam \
    # $RSEM_DIR  \ 
    # $OUTPUT_DIR/RSEM_counts/hcc1395_normal_rep1 >& $OUTPUT_DIR/RSEM_counts/hcc1395_normal_rep1.log # change the file name to match the sample
        
    subprocess.run(["rsem-calculate-expression", "--bam", "--no-bam-output", "-p", "5", "--paired-end", "--forward-prob", "0",
                    os.path.join(alignment_output_dir, name + "Aligned.toTranscriptome.out.bam"), rsem_dir,
                    os.path.join(rsem_output_dir, name), ">&", os.path.join(rsem_output_dir, name + ".log")])


# # Loop through each FASTQ file and run the pipeline
# for fastq_file 
#     sample_name = os.path.splitext(fastq_file)[0]
#     alignment_output_dir = os.path.join(output_directory, sample_name + "_STAR")

#     if not os.path.exists(alignment_output_dir):
#         os.makedirs(alignment_output_dir)
    
#     # Align reads to the reference genome using STAR
#     subprocess.run([star_path, "--genomeDir", reference_genome_dir, 
#                     "--readFilesIn", os.path.join(data_directory, fastq_file),
#                     "--outFileNamePrefix", alignment_output_dir + "/",
#                     "--runThreadN", "NumberOfThreads",
#                     "--outSAMtype", "BAM", "SortedByCoordinate",
#                     "--outSAMunmapped", "Within",
#                     "--outSAMattributes", "Standard"])
    
#     # The output BAM file from STAR is already sorted, so we only need to index it
#     bam_output_path = os.path.join(alignment_output_dir, "Aligned.sortedByCoord.out.bam")
#     subprocess.run([samtools_path, "index", bam_output_path])
    
#     # Assemble transcripts and estimate abundances using StringTie
#     gtf_output_path = os.path.join(output_directory, f"{sample_name}.gtf")
#     subprocess.run([stringtie_path, bam_output_path, "-G", annotation_file_path, "-o", gtf_output_path])

print("RNA sequencing analysis completed with STAR alignment.")