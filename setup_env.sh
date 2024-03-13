#!/bin/bash

# Define the name of the new environment
ENV_NAME="rna_seq_env"

# Path to the requirements.txt file
REQUIREMENTS_PATH="requirements.txt"

# Check if the requirements.txt file exists
if [ ! -f "$REQUIREMENTS_PATH" ]; then
    echo "The requirements.txt file does not exist. Exiting..."
    exit 1
fi

# Create a new conda environment with Python
echo "Creating new Conda environment named $ENV_NAME with Python..."
conda create -n $ENV_NAME python=3.8 -y

# Activate the new environment
echo "Activating the environment..."
source activate $ENV_NAME

# Install bioinformatics tools and any other required packages from requirements.txt
echo "Installing packages from $REQUIREMENTS_PATH..."
conda install -c bioconda --file $REQUIREMENTS_PATH -y

echo "Setup completed. Activate the environment using 'conda activate $ENV_NAME'"