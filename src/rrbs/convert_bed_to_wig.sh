#!/bin/bash

cd /home/tassilo/repos/embryo_epigenetics

# prefixes=$(cat data/prefixes.txt)
gzip -d rrbs_data/*.gz


parallel -j 8 '[[ ! -f {.} ]] && mv {} {.}' ::: rrbs_data/*.txt
# parallel -j 8 'awk '"'"'BEGIN { OFS="\t"; } {if($10=="CpG") print $1, $2, $2+1, $8;}'"'"' {} > {.}_CpG.bg' ::: rrbs_data/*.bed
LC_COLLATE=C

# TODO: check if it is important that the header goes missing in the below line
# parallel -j 8 '[[ ! -f {.}_CpG.bg ]] && awk '"'"'BEGIN { OFS="\t"; } {if($10=="CpG") print $1, $2, $2+1, $8;}'"'"' {} | sort -k1,1 -k2,2n > {.}_CpG.bg' ::: rrbs_data/*.bed
parallel -j 8 '
if [[ ! -f {.}_CpG.bedgraph ]]; then
    HEADER="chromosome\tstart\tend\tmethylation"
    # First create the content without header
    awk '"'"'BEGIN { OFS="\t"; }
        {if($10=="CpG") print $1, $2, $2+1, $8;}'"'"' {} | \
    sort -k1,1 -k2,2n > {.}_CpG.tmp

    # Add header and content to final file
    echo -e "$HEADER" > {.}_CpG.bedgraph
    cat {.}_CpG.tmp >> {.}_CpG.bedgraph
    rm {.}_CpG.tmp
fi' ::: rrbs_data/*.bed
parallel -j 8 '[[ ! -f {.}.bw ]] && bedGraphToBigWig {.}.bedgraph hg19_chrom.sizes {.}.bw' ::: rrbs_data/*_CpG.bedgraph



for cell_type in Sperm Zygote Oocyte "2-cell" "4-cell" "8-cell" Morula ICM TE; do
    for meth_type in "CpG"; do
        # Get all replicates for this condition
        files=$(ls rrbs_data/*${cell_type}*${meth_type}*.bw)
        echo "$files"
        echo

        if [ ${#files[@]} -gt 0 ]; then
            echo "Processing ${cell_type} ${meth_type}"
            echo "run command: wiggletools ratio sum ${files[@]} : sum map unit map offset 1 map default -1 ${files[@]} > rrbs_data/${cell_type}_${meth_type}.wig"
        wiggletools ratio sum "${files[@]}" : sum map unit map offset 1 map default -1 "${files[@]}" > "rrbs_data/${cell_type}_${meth_type}.wig"
        wig2bed < "rrbs_data/${cell_type}_${meth_type}.wig" > "rrbs_data/${cell_type}_${meth_type}.bed"
        fi
    done
done

# for prefix in "${prefixes[@]}"; do
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
