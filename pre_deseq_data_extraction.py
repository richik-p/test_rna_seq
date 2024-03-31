
# generate a file in R to run DESeq2

# Get the TPM_counts column from each os.path.join(output_directory, "RSEM_counts", "${name}") file and write it to a new file


import gen_bash
import os
import pandas as pd
import argparse


output_directory = gen_bash.output_directory

list_of_rsem_files = [f for f in os.listdir(os.path.join(output_directory, "RSEM_counts")) if f.endswith('genes.results')]
list_of_rsem_files.sort()
# Get the TPM_counts column from each os.path.join(output_directory, "RSEM_counts", "${name}") file and write it to a new file

# Create a new text file called TPM_counts.txt to write the TPM counts to
# Write gene_id from one file to the new file
# write all the TPM counts from the other files to the new file after the gene_id as tab separated values

# Get the gene_id column from the first file

first_file = list_of_rsem_files[0]
first_file_path = os.path.join(output_directory, "RSEM_counts", first_file)
first_file_df = pd.read_csv(first_file_path, sep='\t')
gene_id = first_file_df['gene_id']

# get the TPM_counts column from the other files and accumulate them into a dataframe including the gene_id column

TPM_counts_df = pd.DataFrame(gene_id)
raw_counts_df = pd.DataFrame(gene_id)
for file in list_of_rsem_files:
    file_path = os.path.join(output_directory, "RSEM_counts", file)
    df = pd.read_csv(file_path, sep='\t')
    # name the TPM column after the file without the .genes.results extension
    file = file.replace('.genes.results', '')
    TPM_counts_df[file] = df['TPM']
    raw_counts_df[file] = df['expected_count']

TPM_counts_df.to_csv(os.path.join(output_directory, "TPM_counts.txt"), sep='\t', index=False)
raw_counts_df.to_csv(os.path.join(output_directory, "raw_counts.txt"), sep='\t', index=False)