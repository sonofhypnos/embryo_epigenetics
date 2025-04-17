# Set up once in your shell or source file
project_dir="/home/tassilo/embryo_epigenetics"
chrom_sizes="$project_dir/data/hg19.chrom.sizes"

bed_to_bigbed() {
    local input="$1"
    local chrom_sizes="$2"

    if [[ ! -f "$input" || ! -f "$chrom_sizes" ]]; then
        echo "Usage: bed_to_bigbed <input.bed> <chrom.sizes>" >&2
        return 1
    fi

    if ! command -v bedToBigBed &>/dev/null; then
        echo "ERROR: 'bedToBigBed' not found in PATH. Is your conda env activated?" >&2
        return 1
    fi

    local base="${input##*/}"
    local name="${base%%.*}"
    local output="${name}.bb"
    echo "Converting $input to $output"

    local tmp
    tmp=$(mktemp)
    sort -k1,1 -k2,2n "$input" >"$tmp"
    bedToBigBed "$tmp" "$chrom_sizes" "$output"
    rm "$tmp"

    echo "Wrote output to $output"
}

bed_to_bigwig() {
    local input="$1"
    local chrom_sizes="$2"

    if [[ ! -f "$input" || ! -f "$chrom_sizes" ]]; then
        echo "Usage: bed_to_bigwig <input.bedGraph> <chrom.sizes>" >&2
        return 1
    fi

    if ! command -v bedGraphToBigWig &>/dev/null; then
        echo "ERROR: 'bedGraphToBigWig' not found in PATH. Is your conda env activated?" >&2
        return 1
    fi

    local base="${input##*/}"
    local name="${base%%.*}"
    local output="${name}.bw"
    echo "Converting $input to $output"

    local tmp
    tmp=$(mktemp)
    sort -k1,1 -k2,2n "$input" >"$tmp"
    bedGraphToBigWig "$tmp" "$chrom_sizes" "$output"
    rm "$tmp"

    echo "Wrote output to $output"
}

