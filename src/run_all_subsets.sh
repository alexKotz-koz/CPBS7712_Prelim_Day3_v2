#!/bin/bash

# Directory where the subset files are located
dir="data/biosample_data/bat/subsets"

echo $dir
ls

# Iterate over each subset file in the directory
for filepath in "$dir"/*_sliding_window_subset_*.fastq
do
    file=$(basename "$filepath")
    
    echo ""
    echo $file
    echo ""
    # Run main.py in the background
    python3 main.py -k 25 -biosample "$file" &
    # Get the PID of the background process
    pid=$!
    echo "Current PID: $pid"
    # Run a sleep command in the background that kills the main.py process after 10 minutes
    { sleep 600; kill $pid; } &
    # Wait for the main.py process to finish
    wait $pid
    # If the main.py process was killed, its exit status will be 137
    if [ $? -eq 137 ]; then
        echo "Execution of main.py on $file exceeded the time limit and was terminated."
    fi
done