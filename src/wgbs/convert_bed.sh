#!/bin/bash

# Input BED file as first argument
input_bed=$1
output_prefix=${input_bed%.*}

# Filter out lambda chromosome and save to new BED file
grep -v "^lambda" "$input_bed" > "${output_prefix}_filtered.bed"

# Sort the filtered BED file
sort -k1,1 -k2,2n "${output_prefix}_filtered.bed" > "${output_prefix}_sorted.bed"

# Download hg19 chromosome sizes if not present
if [ ! -f "hg19.chrom.sizes" ]; then
    echo "Downloading hg19 chromosome sizes..."
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
fi

# Convert BED to bedGraph
bedtools genomecov -i "${output_prefix}_sorted.bed" -g hg19.chrom.sizes -bg > "${output_prefix}.bedgraph"

# Sort the bedGraph file
sort -k1,1 -k2,2n "${output_prefix}.bedgraph" > "${output_prefix}_sorted.bedgraph"

# Convert to BigWig
bedGraphToBigWig "${output_prefix}_sorted.bedgraph" hg19.chrom.sizes "${output_prefix}.bw"

# Clean up intermediate files
#rm "${output_prefix}_filtered.bed" "${output_prefix}_sorted.bed" "${output_prefix}.bedgraph" "${output_prefix}_sorted.bedgraph"

echo "BigWig file created: ${output_prefix}.bw"
