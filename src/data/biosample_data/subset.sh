#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: ./subset.sh <input_fastq> <output_fastq>"
    exit 1
fi

total_lines=$(wc -l < "$1")

total_reads=$((total_lines/4))

echo "Total lines: $total_lines"
echo "Total reads: $total_reads"



# Get the first 100 lines of the input file and write them to the output file
head -n 100 "$1" > "$2"