#!/usr/bin/env python3

import subprocess
import sys
import os

def compute_depth_for_file(bam_file, window_size, step_size, output_file):
    cmd = ["samtools", "depth", "-a", bam_file]
    depth = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    depths, errors = depth.communicate()

    if depth.returncode != 0:
        print(f"Error with samtools depth: {errors}", file=sys.stderr)
        sys.exit(1)

    chrom = ""
    pos = 0
    window_depths = []

    with open(output_file, 'w') as out:
        for line in depths.split("\n"):
            if not line:
                continue
            
            fields = line.split("\t")
            curr_chrom, curr_pos, curr_depth = fields[0], int(fields[1]), int(fields[2])

            if chrom == "":
                chrom = curr_chrom

            # If moving to a different chrom or passed the window, report and slide
            if curr_chrom != chrom or curr_pos >= pos + window_size:
                if window_depths:
                    avg_depth = sum(window_depths) / len(window_depths)
                    out.write(f"{chrom}\t{pos}\t{pos+window_size-1}\t{avg_depth:.2f}\n")
                    pos += step_size
                    window_depths = window_depths[step_size:]
            
            window_depths.append(curr_depth)

            if curr_chrom != chrom:
                chrom = curr_chrom
                pos = curr_pos

        if window_depths:
            avg_depth = sum(window_depths) / len(window_depths)
            out.write(f"{chrom}\t{pos}\t{pos+window_size-1}\t{avg_depth:.2f}\n")

def batch_process_files_in_directory(input_directory, window_size, step_size, output_directory):
    for root, dirs, files in os.walk(input_directory):
        for file in files:
            if file.endswith(".sorted"):
                bam_file = os.path.join(root, file)
                output_file = os.path.join(output_directory, f"{os.path.splitext(file)[0]}_500readdepth.txt")
                compute_depth_for_file(bam_file, window_size, step_size, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <input_directory> <window_size> <step_size> <output_directory>", file=sys.stderr)
        sys.exit(1)

    input_directory = sys.argv[1]
    window_size = int(sys.argv[2])
    step_size = int(sys.argv[3])
    output_directory = sys.argv[4]

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    batch_process_files_in_directory(input_directory, window_size, step_size, output_directory)