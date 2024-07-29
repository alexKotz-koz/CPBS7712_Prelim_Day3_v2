#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: ./subset.sh <input_fastq> <output_fastq>"
    exit 1
fi

total_lines=$(wc -l < "$1")
total_reads=$((total_lines/4))

lines_per_subset=$((total_lines/10))

# Split the input file into 10 equally sized bins
for i in {1..10}
do
    start_line=$(((i-1)*lines_per_subset+1))
    end_line=$((start_line+lines_per_subset-1))

    # Adjust end_line to the line before the next line that starts with '@', indicating 4 lines per read
    end_line=$(awk -v start="$end_line" '$0 ~ /^@/ && NR > start {print NR-1; exit}' "$1")

    sed -n "${start_line},${end_line}p" "$1" > "${2}_total_subset_${i}.fastq"
done

# Window (number of reads(lines)) and step size
window_size=200
step_size=50

# Create 4 subsets of each of the 10 subsets
for i in {1..10}
do
    for j in {1..4}
    do
        start_line=$(((j-1)*step_size+1))
        end_line=$((start_line+window_size-1))

        # Asure start line is the read id (first char is '@')
        start_line=$(awk -v start="$start_line" '$0 ~ /^@/ && NR >= start {print NR; exit}' "${2}_total_subset_${i}.fastq")

        # Adjust end_line to the line before the next line that starts with '@', indicating 4 lines per read
        end_line=$(awk -v start="$end_line" '$0 ~ /^@/ && NR > start {print NR-1; exit}' "${2}_total_subset_${i}.fastq")

        sed -n "${start_line},${end_line}p" "${2}_total_subset_${i}.fastq" > "${2}_total_subset_${i}_sliding_window_subset_${j}.fastq"
    done
done