import sys
import os
import argparse
from itertools import islice


def count_lines(filename):
    with open(filename) as f:
        return sum(1 for _ in f)


def write_subset(filename, start, end, output_filename):
    with open(filename) as f, open(output_filename, "w") as out:
        for line in islice(f, start, end):
            out.write(line)


def adjust_line(filename, start, direction):
    with open(filename) as f:
        if direction == "forward":
            for i, line in enumerate(f, start=1):
                if i >= start and line.startswith("@"):
                    return i
        elif direction == "backward":
            for i, line in reversed(list(enumerate(f, start=1))):
                if i <= start and line.startswith("@"):
                    return i + 1
    return start


def main(input_file, output_file):
    total_lines = count_lines(input_file)
    lines_per_subset = total_lines // 10
    window_size = 200
    step_size = 50

    for i in range(10):
        start_line = i * lines_per_subset
        end_line = start_line + lines_per_subset
        end_line = adjust_line(input_file, end_line, "backward")
        write_subset(
            input_file, start_line, end_line, f"{output_file}_total_subset_{i+1}.fastq"
        )

    for i in range(10):
        subset_lines = count_lines(f"{output_file}_total_subset_{i+1}.fastq")
        start_line = subset_lines - window_size
        for j in range(4):
            end_line = start_line + window_size
            start_line = adjust_line(
                f"{output_file}_total_subset_{i+1}.fastq", start_line, "forward"
            )
            end_line = adjust_line(
                f"{output_file}_total_subset_{i+1}.fastq", end_line, "backward"
            )
            write_subset(
                f"{output_file}_total_subset_{i+1}.fastq",
                start_line,
                end_line,
                f"{output_file}_total_subset_{i+1}_sliding_window_subset_{j+1}.fastq",
            )
            start_line -= step_size


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a FASTQ file into subsets.")
    parser.add_argument("input_file", help="The input FASTQ file.")
    parser.add_argument("output_file", help="The output FASTQ file.")
    args = parser.parse_args()

    main(args.input_file, args.output_file)
