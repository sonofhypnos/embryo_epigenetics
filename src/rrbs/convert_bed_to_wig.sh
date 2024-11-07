#!/bin/bash

# Configuration variables
DATA_DIR="/home/tassilo/repos/embryo_epigenetics"
# INPUT_DIR="rrbs_data"  # Can be changed to point to different input directories
INPUT_DIR="rrbs_data_test"  # Can be changed to point to different input directories
CELL_TYPES=("PB1" "Sperm" "Zygote" "Oocyte" "2-cell" "4-cell" "8-cell" "Morula" "ICM" "TE")
METH_TYPES=("CpG")

# Change to working directory
cd "${DATA_DIR}"

# Set locale for sorting
export LC_COLLATE=C

# Decompress files if needed
gzip -d "${INPUT_DIR}"/*.gz
parallel -j 8 '[[ ! -f {.} ]] && mv {} {.}' ::: "${INPUT_DIR}"/*.txt

    # # Add header and content to final file
    # echo -e "${HEADER}" > {.}_CpG.bedgraph
    # cat {.}_CpG.tmp >> {.}_CpG.bedgraph

# Process bed files to bedGraph
parallel -j 8 '
if [[ ! -f {.}_CpG.bedgraph ]]; then
    HEADER="chromosome\tstart\tend\tmethylation"
    # First create the content without header
    awk '"'"'BEGIN { OFS="\t"; }
        {if($10=="CpG") print $1, $2, $2+1, $8;}'"'"' {} | \
    sort -k1,1 -k2,2n > {.}_CpG.bedgraph
fi' ::: "${INPUT_DIR}"/*.bed

command -v bedGraphToBigWig || echo  "run 'conda activate epi_env' before running this script."

# Convert bedGraph to bigWig
parallel -j 8 '[[ ! -f {.}.bw ]] && bedGraphToBigWig {.}.bedgraph hg19_chrom.sizes {.}.bw' ::: "${INPUT_DIR}"/*_CpG.bedgraph

# Process each cell type and methylation type combination
for cell_type in "${CELL_TYPES[@]}"; do
    for meth_type in "${METH_TYPES[@]}"; do
        # Get all replicates for this condition
        mapfile -t files < <(ls "${INPUT_DIR}"/GS*"${cell_type}"*"${meth_type}"*.bw)

        if [ ${#files[@]} -gt 0 ]; then
            echo "Processing ${cell_type} ${meth_type}"
            echo "run command: wiggletools ratio sum ${files[*]} : sum map unit map offset 1 map default -1 ${files[*]} > ${INPUT_DIR}/${cell_type}_${meth_type}.wig"

            wiggletools ratio sum "${files[@]}" : sum map unit map offset 1 map default -1 "${files[@]}" > "${INPUT_DIR}/${cell_type}_${meth_type}.wig"
            # wig2bed < "${INPUT_DIR}/${cell_type}_${meth_type}.wig" > "${INPUT_DIR}/${cell_type}_${meth_type}.bed"
            input_file="${DATA_DIR}/${INPUT_DIR}/${cell_type}_${meth_type}.wig"
            output_file="${DATA_DIR}/${INPUT_DIR}/${cell_type}_${meth_type}.bedgraph"
            wig2bed < "${input_file}" > "${output_file}"

            # Convert to bed and sort/merge overlapping regions
            wig2bed < "${input_file}" | \
                sort -k1,1 -k2,2n | \
                bedtools merge -i stdin -c 4 -o mean > "${output_file}"

            # Convert merged bed back to wig
            # awk 'BEGIN{print "fixedStep chrom="$1" start="$2" step=1"}
            #      {for(i=$2; i<$3; i++) print $4}' "${INPUT_DIR}/${cell_type}_${meth_type}.bed" > "${INPUT_DIR}/${cell_type}_${meth_type}_merged.wig"

            # Now convert to bigwig
            bedToBigWig "${INPUT_DIR}/${cell_type}_${meth_type}.bedgraph" hg19_chrom.sizes "${INPUT_DIR}/${cell_type}_${meth_type}_merged.bw"

            # Clean up intermediate files if desired
            # rm "${INPUT_DIR}/${cell_type}_${meth_type}.wig" "${INPUT_DIR}/${cell_type}_${meth_type}_merged.wig"
        else
            echo "No files found for ${cell_type} ${meth_type}"
        fi
    done
done
# done

# for cell_type in Sperm Oocyte Zygote "2cell" "4cell" "8cell" Morula ICM TE hESC "alpha-Amanitin"; do
#     for meth_type in "ACG.TCG" "GCA.GCC.GCT"; do
#         # Get all replicates for this condition
#         files=($(ls *${cell_type}*${meth_type}*.bw))
# 	echo $files

#         if [ ${#files[@]} -gt 0 ]; then
#             echo "Processing ${cell_type} ${meth_type}"
# 	    echo "run command: wiggletools ratio sum ${files[@]} : sum map unit map offset 1 map default -1 ${files[@]} > ${cell_type}_${meth_type}.wig"
# 	    wiggletools ratio sum "${files[@]}" : sum map unit map offset 1 map default -1 "${files[@]}" > "${cell_type}_${meth_type}.wig"
# 	    wig2bed < "${cell_type}_${meth_type}.wig" > "${cell_type}_${meth_type}.bed"
#         fi
#     done
# done
