#!/bin/bash
threads=20
ref="GRCh38"
# Check if the CSV file argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
fi

csv_file=$1

# Check if the file exists
if [ ! -f "$csv_file" ]; then
    echo "File not found!"
    exit 1
fi

# Read the CSV file line by line
line_number=0
while IFS=, read -r col1 col2 col3 col4; do
    if [ "$line_number" -ne 0 ]; then
        echo "Sample: $col1"
        nice AmpliconSuite-pipeline.py \
	-s $col1 \
	-o $col2\
	-t $threads \
	--fastqs $col3 $col4 \
	--ref $ref \
	--run_AA \
	--run_AC
    fi
    line_number=$((line_number + 1))
done < "$csv_file"
