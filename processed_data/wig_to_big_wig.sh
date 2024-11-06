#!/bin/bash

# Check if required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input.wig chrom.sizes output.bw"
    exit 1
fi

INPUT_WIG=$1
CHROM_SIZES=$2
OUTPUT_BW=$3
TEMP_FILTERED="${INPUT_WIG}.filtered"

# Filter out lambda entries using grep and sed
# First, mark sections to remove
grep -B1 -A1000000 "^fixedStep.*chrom=lambda" "$INPUT_WIG" > lambda_sections.txt 2>/dev/null || true
grep -B1 -A1000000 "^variableStep.*chrom=lambda" "$INPUT_WIG" >> lambda_sections.txt 2>/dev/null || true

# Create filtered file excluding lambda sections
if [ -s lambda_sections.txt ]; then
    grep -v -f lambda_sections.txt "$INPUT_WIG" > "$TEMP_FILTERED"
else
    cp "$INPUT_WIG" "$TEMP_FILTERED"
fi

# Remove any empty lines
sed -i '/^$/d' "$TEMP_FILTERED"

# Convert to BigWig
wigToBigWig "$TEMP_FILTERED" "$CHROM_SIZES" "$OUTPUT_BW"

# Check if conversion was successful
if [ $? -eq 0 ]; then
    echo "Successfully converted to BigWig: $OUTPUT_BW"
    # Clean up temporary files
    rm -f "$TEMP_FILTERED" lambda_sections.txt
else
    echo "Error during BigWig conversion"
    exit 1
fi
