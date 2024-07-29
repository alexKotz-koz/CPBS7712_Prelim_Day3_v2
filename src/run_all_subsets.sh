#!/bin/bash
original_dir=$(pwd)
cd data
./rm_logs.sh
cd "$original_dir"

# Directory where the subset files are located, can change biosample_data subdir
dir="data/biosample_data/cryoconite/small_subsets"

# Function to check if a process is running. Use in for loop when a time constraint per run is required
is_running() {
    ps -p $1 > /dev/null 2>&1
}


# Iterate over each subset file in the directory
for filepath in "$dir"/*_sliding_window_subset_*.fastq
do
    file=$(basename "$filepath")
    echo ""
    echo $file
    echo ""
    start_time=$(date +%s)

    # Run main.py in the background, can change size of k
    python3 main.py -k 45 -biosample "$file" &
    # Get the PID of the background process
    pid=$!
    echo "Current PID: $pid"
    
    # Wait for the main.py process to finish
    wait $pid
    
    # uncomment if a time constraint is required
    # If the main.py process was killed, its exit status will be 137
    #if [ $? -eq 137 ]; then
    #    echo "Execution of main.py on $file exceeded the time limit and was terminated." | tee -a src/data/logs/app.log
    #fi
    end_time=$(date +%s)
    execution_time=$((end_time - start_time))

    cd data/logs
    echo "______________________________________________________" | tee -a app.log 
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a app.log
    echo "Total execution time of the script: $execution_time seconds." | tee -a app.log

    # collect results
    cd "$original_dir"
    cd ../  
    cd experimental_results
    filename=$(basename "$file" .fastq)
    if mkdir "$filename"; then
        echo "Directory '$filename' created."
    else
        echo "Failed to create directory '$filename'."
        exit 1
    fi
    cp "$original_dir/data/logs/app.log" "$(pwd)/$filename"
    cp "$original_dir/data/reports/virome_report.txt" "$(pwd)/$filename"
    cp "$original_dir/data/reports/${filename}_tax.csv" "$(pwd)/$filename"
    cp "$original_dir/data/reports/${filename}.pdf" "$(pwd)/$filename"
    cd "$original_dir/data"
    ./rm_logs.sh
    cd "$original_dir"
done

cd ..
cd experimental_results
./collect_results.sh


