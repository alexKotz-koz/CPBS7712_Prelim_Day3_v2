#!/bin/bash
original_dir=$(pwd)
echo $(date '+%Y-%m-%d %H:%M:%S')
cd data
./rm_logs.sh
cd "$original_dir"

# Directory where the subset files are located
dir="data/biosample_data/bat/small_subsets"

echo $dir

# Function to check if a process is running
is_running() {
    ps -p $1 > /dev/null 2>&1
}
start_time=$(date +%s)

# Iterate over each subset file in the directory
for filepath in "$dir"/*_sliding_window_subset_*.fastq
do
    file=$(basename "$filepath")
    
    echo ""
    echo $file
    echo ""
    # Run main.py in the background
    python3 main.py -k 35 -biosample "$file" &
    # Get the PID of the background process
    pid=$!
    echo "Current PID: $pid"
    
    # Run a sleep command in the background that kills the main.py process after 30 minutes
    #(
    #    sleep 1800
    #    if is_running $pid; then
    #        kill $pid
    #    fi
    #) &
    
    # Wait for the main.py process to finish
    wait $pid
    
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

    cd "$original_dir"
    cd ../  # Go up two levels to the root directory
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



