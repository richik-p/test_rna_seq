#!/bin/bash

# add configuration before running with 
# aws configure

# ensure enough memory in the ec2 instance to hold the buckets content

# INSTALLATION
# install python
yes | sudo yum install python # how to automate with y
# install pip
curl -O https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py --user
# install git
sudo yum update -y
sudo yum install git -y
# confirm git installation
if git --version &> /dev/null; then echo "git installed"; else echo "git NOT installed, RE-INSTALL";

# BUCKET SETTINGS
# install aws cli
wget https://s3.amazonaws.com/mountpoint-s3-release/latest/x86_64/mount-s3.rpm
yes | sudo yum install ./mount-s3.rpm

DIR_NAME="test1"
mkdir $DIR_NAME # TODO: set local directory name
aws s3 sync s3://test-bucket-1-r $DIR_NAME # TODO: add bucket name and local directory created

cd $DIR_NAME/

# add environment variable for unzip
export UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE
# unzip the reference genome (takes some time)
unzip hg38_STAR_index.zip # TODO: set the name of reference genome file
mkdir rsem # make directory to contain the zip files
unzip -o rsem_hg38_ref.zip -d rsem # TODO: set the name of reference genome file

# remove the zip files
rm hg38_STAR_index.zip
rm rsem_hg38_ref.zip

# untar the dataset into data
mkdir -p data
tar -xvf RNA_seq_example_set.tar -C data # TODO: change name of the data file

cd ..
source .bashrc

# clone the rna repository
git clone https://github.com/richik-p/test_rna_seq.git

# install conda - takes 8GB memory
wget https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh # check latest version here: https://repo.anaconda.com/archive/
bash Anaconda3-2024.02-1-Linux-x86_64.sh

# get adapter file
mkdir $DIR_NAME/adapter
cd $DIR_NAME/adapter
wget https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa

cd ../..

python test_rna_seq/gen_bash.py --data_dir ./tmp1/data/ --output_dir ./results --ref_genome_dir ./tmp1/hg38_STAR_index --rsem ./tmp1/rsem --th 10 --adap ./tmp1/adapter/TruSeq3-PE-2.fa
