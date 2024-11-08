#!/bin/bash

# Configuration variables
DATA_DIR="/home/tassilo/repos/embryo_epigenetics"
INPUT_DIR="rrbs_data"  # Can be changed to point to different input directories
# INPUT_DIR="rrbs_data_test"  # Can be changed to point to different input directories
CELL_TYPES=("Sperm" "Zygote" "MII_Oocyte" "2-cell" "4-cell" "8-cell" "Morula" "ICM" "TE" "1st_PB" "2nd_PB" "Liver" "PN")
METH_TYPES=("CpG")
SEQ_TYPES=("_RRBS_" "_scRRBS_" "_WGBS_")

# Change to working directory
cd "${DATA_DIR}"

# Set locale for sorting
export LC_COLLATE=C

# Decompress files if needed
gzip -d "${INPUT_DIR}"/*.gz
parallel -j 8 '[[ ! -f {.} ]] && mv {} {.}' ::: "${INPUT_DIR}"/GS*.txt

    # # Add header and content to final file
    # echo -e "${HEADER}" > {.}_CpG.bedgraph
    # cat {.}_CpG.tmp >> {.}_CpG.bedgraph

# Process bed files to bedGraph
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
            mapfile -t files < <(ls "${INPUT_DIR}"/"GS"*"${seq_type}"*"${cell_type}"*"${meth_type}"*.bw)

            if [ ${#files[@]} -gt 0 ]; then
                echo "Processing ${cell_type} ${meth_type}"
                echo "run command: wiggletools ratio sum ${files[*]} : sum map unit map offset 1 map default -1 ${files[*]} > ${INPUT_DIR}/${cell_type}_${meth_type}.wig"

                wiggletools ratio sum "${files[@]}" : sum map unit map offset 1 map default -1 "${files[@]}" > "${INPUT_DIR}/${cell_type}_${meth_type}.wig"
                # wig2bed < "${INPUT_DIR}/${cell_type}_${meth_type}.wig" > "${INPUT_DIR}/${cell_type}_${meth_type}.bed"
                input_file="${DATA_DIR}/${INPUT_DIR}/${cell_type}_${meth_type}.wig"
                bed_file="${DATA_DIR}/${INPUT_DIR}/${cell_type}_${meth_type}.bed"
                bedgraph_file="${DATA_DIR}/${INPUT_DIR}/${cell_type}_${meth_type}.bedgraph"
                wig2bed < "${input_file}" > "${bed_file}"

                # Convert to bed (we remove the index column that I don't know it's origin from )
                sort "${bed_file}" -k1,1 -k2,2n | awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $5}' >  "$bedgraph_file"
                # FIXME: check if below line still needed
                # FIXME: next time learn how we can check inbetween if our files are correctly formatted.
                # bedtools merge -i tmp.bed -c 4 -d 0 -o max > "$bedgraph_file"
                #
                    # bedtools merge -i stdin -c 4 -o mean > "${output_file}"
                # awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $5}' "$bed_file" > tmp.bed

                # Convert merged bed back to wig
                # awk 'BEGIN{print "fixedStep chrom="$1" start="$2" step=1"}
                #      {for(i=$2; i<$3; i++) print $4}' "${INPUT_DIR}/${cell_type}_${meth_type}.bed" > "${INPUT_DIR}/${cell_type}_${meth_type}_merged.wig"

                # Now convert to bigwig
                bedGraphToBigWig "${bedgraph_file}" hg19_chrom.sizes "${INPUT_DIR}/${cell_type}_${meth_type}_merged.bw"
                rm tmp.bed

                # Clean up intermediate files if desired
                # rm "${INPUT_DIR}/${cell_type}_${meth_type}.wig" "${INPUT_DIR}/${cell_type}_${meth_type}_merged.wig"
            else
                echo "No files found for ${cell_type} ${meth_type}"
            fi
        done
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
