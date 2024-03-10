
import os
import subprocess
import argparse
import gzip

# Rest of the argparse setup

args = parser.parse_args()

# Function to unzip .gz files
def unzip_gz_file(gz_path, output_path):
    with gzip.open(gz_path, 'rt') as gz_file:
        with open(output_path, 'w') as out_file:
            for line in gz_file:
                out_file.write(line)

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
        unzip_gz_file(gz_file_path, output_file_path)
        # Replace the .gz file in the list with the uncompressed file
        fastq_files[fastq_files.index(gz_file)] = gz_file[:-3]
else:
    # Original listing for non-compressed files
    fastq_files = [f for f in os.listdir(data_directory) if f.endswith('.fastq') or f.endswith('.fq')]

# Continue with the original script for alignment and analysis

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description='RNA-Seq analysis pipeline using STAR for alignment.')
parser.add_argument('--reference_genome_dir', required=True, help='Directory of the STAR reference genome indices')
parser.add_argument('--annotation_file_path', required=True, help='Path to the annotation GTF file')

# Parse the command-line arguments
args = parser.parse_args()

# Assign arguments to variables
reference_genome_dir = args.reference_genome_dir
annotation_file_path = args.annotation_file_path

# Define paths to your tools
star_path = "STAR"
stringtie_path = "stringtie"
samtools_path = "samtools"

# Define the directory containing your FASTQ files
data_directory = "/path/to/your/data/folder"

# Define the output directory for your results
output_directory = "/path/to/your/output/folder"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# List all FASTQ files in the data directory
fastq_files = [f for f in os.listdir(data_directory) if f.endswith('.fastq') or f.endswith('.fq')]

# Loop through each FASTQ file and run the pipeline
for fastq_file in fastq_files:
    sample_name = os.path.splitext(fastq_file)[0]
    alignment_output_dir = os.path.join(output_directory, sample_name + "_STAR")

    if not os.path.exists(alignment_output_dir):
        os.makedirs(alignment_output_dir)
    
    # Align reads to the reference genome using STAR
    subprocess.run([star_path, "--genomeDir", reference_genome_dir, 
                    "--readFilesIn", os.path.join(data_directory, fastq_file),
                    "--outFileNamePrefix", alignment_output_dir + "/",
                    "--runThreadN", "NumberOfThreads",
                    "--outSAMtype", "BAM", "SortedByCoordinate",
                    "--outSAMunmapped", "Within",
                    "--outSAMattributes", "Standard"])
    
    # The output BAM file from STAR is already sorted, so we only need to index it
    bam_output_path = os.path.join(alignment_output_dir, "Aligned.sortedByCoord.out.bam")
    subprocess.run([samtools_path, "index", bam_output_path])
    
    # Assemble transcripts and estimate abundances using StringTie
    gtf_output_path = os.path.join(output_directory, f"{sample_name}.gtf")
    subprocess.run([stringtie_path, bam_output_path, "-G", annotation_file_path, "-o", gtf_output_path])

print("RNA sequencing analysis completed with STAR alignment.")
