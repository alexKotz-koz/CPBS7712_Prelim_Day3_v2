#!/bin/bash

# Define the directories
dir1="logs"
dir2="output_data"

# Delete all files in directory 1
rm -f "$dir1"/*

# Delete all files in directory 2
rm -f "$dir2"/*