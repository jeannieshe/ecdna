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
while IFS=, read input_dir sample output_dir; do
    echo "--------------------------------------------------"
    echo "Input dir: $input_dir, Sample name: $sample, Output dir: $output_dir"
    if [ "$input_dir" != "input-dir" ]; then  # Skip header line explicitly
        echo "Sample: $sample"
        echo "Running script: /home/jeannie/group/common_scripts/AmpliconArchitect/concat_fastqs.sh \"$input_dir\" \"$sample\" -o \"$output_dir\""
        /home/jeannie/group/common_scripts/AmpliconArchitect/concat_fastqs.sh "$input_dir" "$sample" -o "$output_dir"
    fi
done < "$csv_file"

echo "Finished processing file: $csv_file"
