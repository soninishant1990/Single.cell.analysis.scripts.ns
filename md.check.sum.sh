#!/bin/bash

# Output CSV file
output_file="checksums.csv"

# Header for CSV file
echo "Filename,MD5Checksum" > "$output_file"

# Iterate over folders
for folder in /broad/hptmp/hptmp_soni/GBM.FILE/*/; do
    # Iterate over .fastq.gz files in the folder
    for file in "$folder"*.fastq.gz; do
        # Calculate MD5 checksum
        checksum=$(md5sum "$file" | awk '{ print $1 }')
        # Get filename without path
        filename=$(basename "$file")
        # Append to CSV file
        echo "$filename,$checksum" >> "$output_file"
    done
done