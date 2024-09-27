# Create loompy files from fastq files
# Original pipeline: https://linnarssonlab.org/loompy/kallisto/index.html#:~:text=This%20section%20introduces%20the%20loompy,files%20directly%20from%20fastq%20files

# Set Up:
# Install kallisto
# https://pachterlab.github.io/kallisto/
# https://pachterlab.github.io/kallisto/source

# Make tab delimted metadata file with name, technology, and targetnumcells
# Must have metadata.tab for code to know what fastq files to use and such
# Example:
# name    technology  targetnumcells
# 1kPBMC  10xv3       1000

# Set up variables
org = "mouse"

# Where this RNA Velocity Code is:
projDir = 
metadata = 

# Where you're big files are (proj_dir holds the index and packages, fastqs holds the fast files you will run this on):
proj_dir = 
fastqs = 

# Download an index for human data (or build your own for mouse data)
if org == "human":
    # Download an index for human data
    index = proj_dir + "human_GRCh38_gencode.v31.600/"
else:
    # Or Build your own for mouse
    # https://github.com/linnarsson-lab/loompy/tree/6f87b6a3b6fcc7dcefc92532b93369717479b2cd/kallisto
    # This is actually annoying difficult to build your own, I suggest referencing my files:
    index = '/broad/hptmp/ughetta/Komal/'


import csv
sampleIDs = []

## ==================== ##
##  Using module 'csv'  ##
## ==================== ##
with open(metadata) as to_read:
    reader = csv.reader(to_read, delimiter = "\t")
    desired_column = [0] # SampleID column

    for row in reader:     # read one row at a time
        sampleIDs += list(row[i] for i in desired_column)   # build the output row (process)
                         
sampleIDs = sampleIDs[1:] # Skips header
sampleIDs

# Run the loompy fromfq command for each sample
# loompy fromfq sampleID.loom sampleID index metadata sampleID1_S1_L001_R1_001.fastq.gz sampleID1_S1_L001_R2_001.fastq.gz sampleID1_S1_L002_R1_001.fastq.gz sampleID1_S1_L002_R2_001.fastq.gz

# Generate shell sript for job submission to run loompy for each sample
outfile = projDir + "fastq2loom_job.sh"
f = open(outfile, "w")

f.write("#!/bin/bash\n")
f.write("#$ -M ughetta@login.broadinstitute.org\n")
f.write("#$ -N job.fastq2loom\n")
f.write("#$ -cwd\n")
f.write("#$ -q broad\n")
f.write("#$ -l h_vmem=8g\n")
f.write("#$ -pe smp 1\n")
f.write("#$ -binding linear:1\n")
f.write("#$ -l h_rt=12:00:00\n")
f.write("#$ -e fastq2loom.err\n\n")
f.write("#$ -M ughetta@broadinstitute.org\n")
f.write("#$ -m bea\n")

f.write("source /broad/software/scripts/useuse\n")
#f.write("use .anaconda3-5.3.1\n")
f.write("use UGER\n")
f.write("use Anaconda3\n")
f.write("use R-4.1\n")
f.write("source activate /home/unix/ughetta/conda/scrnatools\n")
f.write("use UGER\n")
f.write("use Anaconda3\n")
f.write("use R-4.1\n\n")

#fixes libhdf5.so.10 error
f.write("export PATH=$HOME/bin:$PATH\n\n")
f.write("export LD_LIBRARY_PATH=$HOME/lib/:$LD_LIBRARY_PATH\n\n")

f.write("sampleID=${1}\n")
f.write("index=${2}\n")
f.write("metadata=${3}\n")
f.write("fastqs=${4}\n")
f.write("out_dir=${5}\n")

f.write("echo $sampleID\n")
f.write("echo $index\n")
f.write("echo $metadata\n")
f.write("echo $fastqs\n")
f.write("echo $out_dir\n")

f.write("cd /broad/hptmp/ughetta/Komal/\n")
f.write("pwd\n")

# EDIT THIS LINE TO MATCH FASTQ NAMING: this is an example of what I used
f.write("loompy fromfq ${fastqs}${sampleID}.loom ${sampleID} ${index} ${metadata} ${fastqs}KOAL01_Centa${sampleID}_0_G_S*_L001_R1_001.fastq.gz ${fastqs}KOAL01_Centa${sampleID}_0_G_S*_L001_R2_001.fastq.gz ${fastqs}KOAL01_Centa${sampleID}_0_G_S*_L002_R1_001.fastq.gz ${fastqs}KOAL01_Centa${sampleID}_0_G_S*_L002_R2_001.fastq.gz\n")

f.write("echo done\n")
f.close()

# Make sh file executable
import os
os.system("chmod +x " + outfile)

# Checking premissions for sanity check
# os.system("ls -l " + outfile)

import os.path
from os import path

## Submit job on server for files that do not exist
for ID in sampleIDs:
    if path.exists(fastqs+ID+".loom") == False:
        os.system('qsub ' + projDir + 'fastq2loom_job.sh ' + ID + ' ' + index + ' ' + metadata + ' '+ fastqs + ' ' + proj_dir)      