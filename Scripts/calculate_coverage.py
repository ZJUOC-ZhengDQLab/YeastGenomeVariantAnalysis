import os
import subprocess

def calculate_depth(file_path, condition):
    """
    Calculate the average coverage under specific conditions using the samtools depth command.
    :param file_path: Path to the input BAM file
    :param condition: awk condition to filter regions of interest
    :return: Calculated average coverage (float)
    """
    try:
        # Construct the samtools and awk command
        command = f"samtools depth {file_path} | awk '{condition}'"
        
        # Use subprocess to get the command output
        result = subprocess.check_output(command, shell=True, text=True)
        
        # Convert the result to a float
        return float(result.strip()) if result.strip() else 0.0
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        return 0.0

def main():
    """
    Main function: Iterate through files ending with .sorted in the current directory and calculate coverage.
    The results are saved in the current directory, with each sample generating a corresponding result file.
    """
    # Get all files ending with .sorted in the current directory
    current_dir = os.getcwd()
    sorted_files = [f for f in os.listdir(current_dir) if f.endswith(".sorted")]

    if not sorted_files:
        print("No .sorted files found in the current directory.")
        return

    for sorted_file in sorted_files:
        # Sample name is the filename without the extension
        sample_name = os.path.splitext(sorted_file)[0]

        # Define the output filename
        output_file = os.path.join(current_dir, f"{sample_name}_coverage.txt")

        # Calculate chromosome coverage (excluding chr=mt and chr=2μ)
        chromosome_condition = '$1 != "chr=mt" && $1 != "chr=2micro" && $3 > 0 {sum+=$3; count++} END {if (count > 0) print sum/count; else print 0}'
        chromosome_depth = calculate_depth(sorted_file, chromosome_condition)

        # Calculate mitochondrial coverage
        mt_condition = '$1 == "chr=mt" && $3 > 0 {sum+=$3; count++} END {if (count > 0) print sum/count; else print 0}'
        mt_depth = calculate_depth(sorted_file, mt_condition)

        # Calculate 2μ chromosome coverage
        twomicro_condition = '$1 == "chr=2micro" && $3 > 0 {sum+=$3; count++} END {if (count > 0) print sum/count; else print 0}'
        twomicro_depth = calculate_depth(sorted_file, twomicro_condition)

        # Write the results to the file (no header)
        with open(output_file, "w") as f:
            f.write(f"{sample_name}\t{chromosome_depth}\t{mt_depth}\t{twomicro_depth}\n")

        # Notify the user that the file was successfully generated
        print(f"Coverage results saved to {output_file}")

if __name__ == "__main__":
    main()