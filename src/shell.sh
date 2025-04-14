# Functions to import into shell if needed
project_dir="/home/tassilo/embryo_epigenetics"
chrom_sizes="$project_dir/data/hg19.chrom.sizes"

bed_to_bigbed() {
    local input="$1"
    local chrom_sizes="$2"

    if [[ ! -f "$input" || ! -f "$chrom_sizes" ]]; then
        echo "Usage: bedgraph_to_bigbed <input.bedGraph|.bed> <chrom.sizes>" >&2
        return 1
    fi

    local base="${input##*/}" # Strip path
    local name="${base%%.*}"  # Strip extension
    local output="${name}.bb"
    echo "Converting $input to $output"

    grep -v '^track' "$input" | sort -k1,1 -k2,2n >temp.bed
    bedToBigBed temp.bed "$chrom_sizes" "$output"
    rm temp.bed

    echo "Wrote output to $output"
}

bed_to_bigwig() {
    local input="$1"
    local chrom_sizes="$2"

    if [[ ! -f "$input" || ! -f "$chrom_sizes" ]]; then
        echo "Usage: bed_to_bigwig <input.bedGraph|.bed> <chrom.sizes>" >&2
        return 1
    fi

    local base="${input##*/}" # Strip path
    local name="${base%%.*}"  # Strip extension
    local output="${name}.bw"
    echo "Converting $input to $output"

    grep -v '^track' "$input" | sort -k1,1 -k2,2n >temp.bed
    bedGraphToBigWig temp.bed "$chrom_sizes" "$output"
    rm temp.bed

    echo "Wrote output to $output"
}
