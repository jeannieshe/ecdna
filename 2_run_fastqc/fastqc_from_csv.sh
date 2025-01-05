#!/bin/bash

# Check if the CSV file argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
fi

csv_file="$1"

# Check if the file exists
if [ ! -f "$csv_file" ]; then
    echo "File not found: $csv_file"
    exit 1
fi

# Read the CSV file line by line, skipping the header
# Make sure the CSV file ends in a new line character!
while IFS=, read input_dir; do
    echo "--------------------------------------------------"
    if [ "$input_dir" != "input-dir" ]; then  # Skip header line explicitly
        echo "Input dir: $input_dir"
        file1=($(find "$input_dir" -type f -name '*1_merged.fq.gz') )
        file2=($(find "$input_dir" -type f -name '*2_merged.fq.gz') )

        echo "Running script: fastqc \"$file1\" \"$file2\""
        fastqc "$file1" "$file2"
    fi
done < "$csv_file"

echo "Finished processing file: $csv_file"
