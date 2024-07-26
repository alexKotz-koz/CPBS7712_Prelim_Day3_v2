if grep -q '^$' biosample_data/biofilm/SRR24581281.fastq; then
    echo "There are empty lines in the file."
else
    echo "There are no empty lines in the file."
fi