#!/bin/bash

FILE="deseq.yml"
# Extract the name of the environment from the FILE file
ENV_NAME=$(grep 'name:' $FILE | cut -d ' ' -f 2)

# Check if the environment already exists
conda info --envs | grep "^$ENV_NAME" > /dev/null

if [ $? -ne 0 ]; then
    echo "Environment '$ENV_NAME' does not exist. Creating..."
    conda env create -f $FILE
    echo "Environment '$ENV_NAME' created."
else
    echo "Environment '$ENV_NAME' already exists."
fi

# Activate the environment
echo "Activating environment '$ENV_NAME'."
conda init
conda activate $ENV_NAME