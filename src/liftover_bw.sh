#!/bin/bash -
#title          :liftover_bw.sh
#description    :Liftover big wig file and transform back to bigwig file
#author         :Tassilo Neubauer
#date           :20250320
#version        :0.1
#usage          :./liftover_bw.sh
#notes          :
#bash_version   :5.1.16(1)-release
#============================================================================

#!/bin/bash

# Check if a filename was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <wigfile>"
    exit 1
fi

# Get the input file name
WIGFILE="$1"

# Check if file has an extension
if [[ "$WIGFILE" != *.* ]]; then
    echo "Error: Input file must have an extension (e.g., .wig, .bedgraph, etc.)"
    exit 1
fi

# Extract base name without extension for output naming
BASENAME=$(basename "$WIGFILE" | sed 's/\.[^.]*$//')

# Convert WIG to bedGraph, sort, and perform liftOver
wiggletools write_bg - "$WIGFILE" |
    sort -k1,1 -k2,2n >temp.bedGraph &&
    liftOver temp.bedGraph ../hg19ToHg38.over.chain temp_lifted.bedGraph "${BASENAME}_unlifted.txt" &&
    sort -k1,1 -k2,2n temp_lifted.bedGraph >temp_lifted_sorted.bedGraph || echo "Something went wrong in the conversion"

# Merge overlapping regions using bedtools
bedtools merge -i temp_lifted_sorted.bedGraph -c 4 -o mean >temp_lifted_merged.bedGraph

# Convert merged bedGraph to bigWig
bedGraphToBigWig temp_lifted_merged.bedGraph ../hg38.chrom.sizes "${BASENAME}_hg38.bw"

# Clean up temporary files
rm temp.bedGraph temp_lifted.bedGraph temp_lifted_sorted.bedGraph temp_lifted_merged.bedGraph

echo "Conversion complete: ${BASENAME}_hg38.bw"
