#!/bin/bash
# Removes all log files
# Define the directories
dir1="logs"
dir2="output_data"
dir3="reports"

# Delete all files in directory 1
rm -f "$dir1"/*

# Delete all files in directory 2
rm -f "$dir2"/*

rm -f "$dir3"/*