#!/bin/bash

# Get the current directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Create a temporary file to store all the data
temp_file=$(mktemp)

# Create a new file to store the results
result_file="$DIR/result.csv"

# Create a new file to store the execution times
execution_stats_file="$DIR/execution_stats.txt"

# Variable to handle header
header_handled=false

# Variable to store total execution time
total_execution_time=0

# Variable to store the number of reads in the original biosample file
original_biosample_reads=0

# Iterate through each subdirectory
for subdir in "$DIR"/*; do
    if [ -d "$subdir" ]; then
        cd "$subdir"
        # Iterate over each CSV file in the subdirectory
        for csv_file in *_tax.csv; do
            if [ -f "$csv_file" ]; then
                # If header is not handled, copy the header from the first file and mark it as handled
                if ! $header_handled ; then
                    head -n 1 "$csv_file" > "$result_file"
                    header_handled=true
                fi
                # Append the data (excluding the header) from the CSV file to the temporary file
                tail -n +2 "$csv_file" >> "$temp_file"
            fi
        done
        # Extract the execution time from the app.log file and add it to the total
        if [ -f app.log ]; then
            execution_time=$(awk '/Total execution time of the script:/ {print $7}' app.log)
            total_execution_time=$((total_execution_time + execution_time))
        fi
        # Extract the number of reads in the original biosample file from the virome_report.txt file
        if [ -f virome_report.txt ] && [ $original_biosample_reads -eq 0 ]; then
            original_biosample_reads=$(awk -F': ' '/Number of Reads in Original Biosample File/ {print $2}' virome_report.txt)
        fi
    fi
done

# Sort the data by virus name, sum the abundance for each virus, and append the results to the result file
awk -F, '{a[$1] += $15; lines[$1]=$0} END{for (i in a) {split(lines[i],fields,","); fields[15]=a[i]; $0=fields[1]; for(j=2;j<=length(fields);j++){$0=$0","fields[j]}; print $0}}' "$temp_file" >> "$result_file"

# Write the total execution time to the execution_stats.txt file
echo "Total execution time: $total_execution_time seconds" > "$execution_stats_file"

# Calculate the percentage of reads used and write it to the execution_stats.txt file
percentage_of_reads_used=$(echo "scale=10; (40/$original_biosample_reads)*100" | bc)
echo "Percentage of reads used: $percentage_of_reads_used%" >> "$execution_stats_file"

# Remove the temporary file
rm "$temp_file"