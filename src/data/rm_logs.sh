#!/bin/bash
# Removes all log files

logs="logs"
output_data="output_data"
reports="reports"

# Delete all files in directory 1
rm -f "$logs"/*

# Delete all files in directory 2
rm -f "$output_data"/*

rm -f "$reports"/*