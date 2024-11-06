#!/bin/bash

cd ../../
conda activate epi_env

# prefixes=$(cat data/prefixes.txt)
gzip rrbs_data/*.gz


parallel -j 8 'mv {} {.}' ::: rrbs_data/*.txt
parallel -j 8 'awk '"'"'BEGIN { OFS="\t"; } {if($10=="CpG") print $1, $2, $2+1, $8;}'"'"' {} > {.}_CpG.bg' ::: rrbs_data/*.bed
parallel -j 8 'bedGraphToBigWig {.}.bedgraph hg19_chrom.sizes {.}.bw' ::: rrbs_data/*_CpG.bg


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
