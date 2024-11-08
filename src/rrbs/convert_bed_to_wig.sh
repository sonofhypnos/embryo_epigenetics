#!/bin/bash

DATA_DIR="/home/tassilo/repos/embryo_epigenetics"
INPUT_DIR="$DATA_DIR/rrbs_data"  # Default input directory

# Change directory to test directory based on arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --test)
            INPUT_DIR="${INPUT_DIR}_test"
            shift
            ;;
        *)
            shift
            ;;
    esac
done

CELL_TYPES=("Sperm" "Zygote" "MII_Oocyte" "2-cell" "4-cell" "8-cell" "Morula" "ICM" "TE" "1st_PB" "2nd_PB" "Liver" "PN")
METH_TYPES=("CpG")
SEQ_TYPES=("_RRBS" "_scRRBS" "_WGBS")

# Change to working directory
cd "${DATA_DIR}"

# Set locale for sorting
export LC_COLLATE=C

# Decompress files if needed
gzip -d "${INPUT_DIR}"/*.gz
parallel -j 8 '[[ ! -f {.} ]] && [[ -f {} ]] && mv {} {.}' ::: "${INPUT_DIR}"/GS*.txt
parallel -j 8 '[[ ! -f {.} ]] && [[ -f {} ]] && mv {} {.}' ::: "${INPUT_DIR}"/GS*.txt

# Process bed files to bedGraph
# NOTE: in order to merge adjacent windows with bedtools merge, unexplicably the correct value is 1?
parallel -j 8 '
if [[ ! -f {.}_CpG.bedgraph ]]; then
    echo processing {} ...
    HEADER="chromosome\tstart\tend\tmethylation"
    # First create the content without header
    awk '"'"'BEGIN { OFS="\t"; }
        {if($10=="CpG") print $1, $2, $2+1, $8;}'"'"' {} | \
    sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -d 1 -o max > {.}_CpG.bedgraph
fi' :::    "${INPUT_DIR}"/GS*.bed # NOTE: we are adding GS in the beginning, so we don't match the files we are creating later



command -v bedGraphToBigWig || echo  "run 'conda activate epi_env' before running this script."

# Convert bedGraph to bigWig
parallel -j 8 '[[ ! -f {.}.bw ]] && echo transforming {} to bigwig file && bedGraphToBigWig {.}.bedgraph hg19_chrom.sizes {.}.bw' ::: "${INPUT_DIR}"/GS*_CpG.bedgraph

# Process each cell type and methylation type combination
for cell_type in "${CELL_TYPES[@]}"; do
    for seq_type in "${SEQ_TYPES[@]}"; do
        for meth_type in "${METH_TYPES[@]}"; do
            # Get all replicates for this condition
            mapfile -t files < <(ls "${INPUT_DIR}"/"GS"*"${seq_type}"*"${cell_type}"*"${meth_type}"*.bw 2>/dev/null)

            if [ ${#files[@]} -gt 0 ]; then
                echo "Processing ${cell_type} ${meth_type}"
                input_file="${INPUT_DIR}/${cell_type}_${meth_type}${seq_type}.wig"
                bed_file="${INPUT_DIR}/${cell_type}_${meth_type}${seq_type}.bed"
                bedgraph_file="${INPUT_DIR}/${cell_type}_${meth_type}${seq_type}.bedgraph"
                bw_file="${INPUT_DIR}/${cell_type}_${meth_type}${seq_type}_merged.bw"
                echo "run command: wiggletools ratio sum ${files[*]} : sum map unit map offset 1 map default -1 ${files[*]} > ${input_file}"

                wiggletools ratio sum "${files[@]}" : sum map unit map offset 1 map default -1 "${files[@]}" > "${input_file}"
                # wig2bed < "${INPUT_DIR}/${cell_type}_${meth_type}.wig" > "${INPUT_DIR}/${cell_type}_${meth_type}.bed"
                wig2bed < "${input_file}" > "${bed_file}"

                # Convert to bed (we remove the index column that I don't know it's origin from )
                sort "${bed_file}" -k1,1 -k2,2n | awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $5}' >  "$bedgraph_file"

                # Now convert to bigwig
                bedGraphToBigWig "${bedgraph_file}" hg19_chrom.sizes "${bw_file}"

                # NOTE: this line is just here to check everything is working as expected
                wiggletools "$bw_file" 1>/dev/null || echo "Big wig files are not correct!" >&2
            fi
        done
    done
done
