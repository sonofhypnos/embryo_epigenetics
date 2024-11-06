#!/bin/bash

# Check if input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input.bw>"
    exit 1
fi

# Input bigWig file
input_bw=$1

# Remove the .bw extension from the input filename
basename=$(basename "$input_bw" .bw)
output_bw="${basename}_lambda_only.bw"

# Create unique temporary directory using process ID
tmp_dir=$(mktemp -d -t lambda_filter.XXXXXX)

# Use trap to ensure cleanup on script exit, including if interrupted
trap 'rm -rf "$tmp_dir"' EXIT

# Use process-specific temporary files
temp_bedGraph="${tmp_dir}/temp_${$}.bedGraph"
lambda_bedGraph="${tmp_dir}/lambda_${$}.bedGraph"
chrom_sizes="${tmp_dir}/lambda_${$}.chrom.sizes"

# Extract just the lambda chromosome using bigWigToBedGraph and bedGraphToBigWig
bigWigToBedGraph "$input_bw" "$temp_bedGraph"
grep "^lambda" "$temp_bedGraph" > "$lambda_bedGraph"

# Create chromosome sizes file for lambda
echo -e "lambda\t48502" > "$chrom_sizes"

# Convert back to bigWig
bedGraphToBigWig "$lambda_bedGraph" "$chrom_sizes" "$output_bw"

# Cleanup is handled automatically by the trap
