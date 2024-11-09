#!/bin/bash


# Configuration variables
DATA_DIR="/home/tassilo/repos/embryo_epigenetics"
INPUT_DIR="$DATA_DIR/data"  # Default input directory

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

# Process each cell type and methylation type separately
for cell_type in Sperm Oocyte Zygote "2cell" "4cell" "8cell" Morula ICM TE hESC "alpha-Amanitin"; do
    for meth_type in "ACG.TCG" "GCA.GCC.GCT"; do
        # Get all replicates for this condition
        files=("$(ls ${INPUT_DIR}/*${cell_type}*${meth_type}*.bw)")
        echo "$files"
        
        if [ ${#files[@]} -gt 0 ]; then
            echo "Processing ${cell_type} ${meth_type}"
	    echo "run command: wiggletools ratio sum ${files[@]} : sum map unit map offset 1 map default -1 ${files[@]} > ${INPUT_DIR}/${cell_type}_${meth_type}.wig"
	    wiggletools ratio sum "${files[@]}" : sum map unit map offset 1 map default -1 "${files[@]}" > "${INPUT_DIR}/${cell_type}_${meth_type}.wig"
	    wig2bed < "${INPUT_DIR}/${cell_type}_${meth_type}.wig" > "${INPUT_DIR}/${cell_type}_${meth_type}.bed"
        fi
    done
done

# Generate UCSC hub.txt
# (reuse our previous Python code for this part since it's just text generation)
