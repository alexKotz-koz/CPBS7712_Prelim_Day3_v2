# CPBS Prelims Day 3
> Metagenomic Sequence Assembler and Virome Characterizer 

**Goal:** Develop a general method for characterizing the virome (viral communities) in a metagenomic sample. The input to your method is a short-read fastq file containing sequences from a metagenomic sample. The output to your method should be a report characterizing the virome, which would be useful to researchers, clinicians, and/or public health officials working in pandemic preparedness.

The outline of this project is as follows:
1. Read in the metagenomic sample (.fastq) and known viruses (.fasta) 
    - importBioSample.py
    - importVirus.py
2. Perform Quality Control on the metagenomic sample:
    - Convert existing quality control ascii values to a Phred score and perform pruning of low quality reads
    - qc.py
3. Break the metagenomic sample reads into k-mers (user defined)
    - readsToKmers.py
4. Create a de Bruijn Graph from the pool of metagenomic sample k-mers
    - deBruijnGraph.py
5. Create contigs via a Depth-First Search traverse across the de Bruijn Graph
    - createContigs.py
6. Search the contigs for substrings of the viral sequences
    - searchForViruses.py

**Note on Input files:** Expected file types of input files are .fastq. Each read should contain 4 lines of information:

- Line 1: Read id, expected id's should start with @, such as: @SRR24581281.1 
- Line 2: Read sequence: ACTGGATCTTCAG
- Line 3: Expected format: +SRR24581281.1 1 length=302
- Line 4: Quality Information, should be ascii characters, such as: AAAAA#EEEEEEE

## Installation

OS X & Linux:
1. Clone or download the repository.
2. Set up miniconda environment:
    - If miniconda is not installed on the local machine, please follow the steps outlined here before continuing: [Miniconda installation](https://docs.anaconda.com/free/miniconda/)
    - Once miniconda is installed, create the conda environment by copying this command into a shell (terminal) with an active base conda environment:
        ```sh
        conda env create -f conda_env.yml
        ```
    - Then activate the new conda environment:
        ```sh
        conda activate conda_env
        ```

## Usage example
1. Open a terminal and navigate to /src/:
```sh
cd src
```
2. Use the following command to run the project: 
```sh
python3 main.py [options]
```
Options:
- `-h, --help`: Show help menu
- `-biosample`: Metagenomic biosample file (required)
- `-k`: User defined size of the k-mers (required)

<br>
Example:

1. After navigating to src (cd src), run the following command:
```sh
python3 main.py -biosample synthetic -k 5
```
2. Let the program execute, you will see time stamps for each component. Such as:
```sh
Time Stamp: Reads to Kmers finished in 0.22631192207336426
Time Stamp: DeBruijn Graph finished in 0.1376938819885254
Time Stamp: Create Contigs finished in 22.05834698677063
```
3. Results will be located in src/data/reports



## Virome Report
To characterize the virome present in a given metagenomic sample, several metrics are calculated:
1. Taxonomy of each virus tested
- Located in the csv file under src/data/reports
2. Relative Abundance of the virus in the metagenomic sample.
- Located in the csv file under src/data/reports
3. Number of metagenomic contigs found in the virus.
- Located in the csv file under src/data/reports
4. Quality control metrics: Information regarding the quality control of each metagenomic sample.
- Located in the pdf file under src/data/reports
5. Number of reads in the original biosample file. This is used in combination with item 6 in this list, to demonstrat the total number of reads used in the experiement.
- Located in the pdf file under src/data/reports
6. Experimental Results: Results for each experiment run can be found under experiemental_results. Each 'experiment' is the result of executing run_all_subsets.sh (located under src/) for each of the metagenomic samples used in the development of this pipeline.
- execution_stats.txt: Contains the total execution time of the experiment and the total percentage of reads used, from the original metagenomic sample.
- result.csv: Contains the taxonomic information of all viruses tested against the metagenomic sample and the aggregated relative abundance of each virus in all of the subsets of the original metagenomic sample

### Example of the 2 Virome Report Files:

- Biosample Information: Found in src/data/reports/synthetic_biosample.pdf
- Taxonomy of Viruses Tested: Found in src/data/reports/synthetic_biosample_tax.csv

## Requirements
- Python 3.9 or higher. 
    - This project was developed in python v3.11 and has been tested with python 3.9 thru 3.12.
- Miniconda (see Installation section for further instructions).
- macOS or Linux based operating system.

## Automation and Supplemental Files:
- Number of total reads in the original biosample files were calculated by:
```sh
grep -o '@' filname | wc -l
```
- Automation Scripts Include:
    - **run_all_subsets.sh**: Executes the metagenomic pipeline on all 40 subsets of the metagenomic sample, generating results for each run in the experiemental_results directory at the root of the project.
        - Note on the instructions below: Please change to the script directory (cd ...) before executing these scripts.
        - To execute this script, first follow the download sequence data instructions (located in the DownloadSequenceData_Instructions.md file). Then adjust the file path for the biosample directory on line 8 of the run_all_subsets.sh. 
        - Then make the run_all_subsets.sh an executable with global permissions:
        ```sh
        chmod +x run_all_subsets.sh
        ```
        - Finally, execute the script:
        ```sh
        ./run_all_subsets.sh
        ```
        - Note this script will take 30 minutes to 2 hours to run, based on the size of the subsets, size of k, and hardware.

    - **subset_biosample.sh**: Creates 40 subsets of the original metagenomic sample via a sliding window technique.
        - To execute this script and create subsets of the original metagenomic sample:
        1. Make the file executable:
        ```sh
        chmod +x subset_biosample.sh
        ```
        2. Create the directories to dump the subsets into:
        - Recommended file structure:
            - src/data/biosample_data/<name of biosample 1>
                - src/data/biosample_data/<name of biosample 1>/small_subsets
        3. Execute the script:
        ```sh
        ./subset_biosample.sh <file location of original biosample.fasta file> <file location of output files>
        ```
        - Example output file arg: src/data/biosample_data/bat/bat1 (bat1 being the prefix to all subset files)
    - **collect_results.sh**: Aggregates the relative abundance for each virus from all 40 subsets and calculates the percentage of reads used in the 40 subsets from the original metagenomic sample.
        - Note: This file assumes that 40 subset fastq files were created, if this number is different, please change the numerator on line 58.
        1. Make the file executable with global permissions:
        ```sh
        chmod +x collect_results.sh
        ```
        2. Execute the file:
        ```sh
        ./collect_results.sh
        ```