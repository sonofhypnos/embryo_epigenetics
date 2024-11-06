#!/bin/bash

cd /home/tassilo/repos/embryo_epigenetics

# prefixes=$(cat data/prefixes.txt)
gzip -d rrbs_data/*.gz


parallel -j 8 '[[ ! -f {.} ]] && mv {} {.}' ::: rrbs_data/*.txt
# parallel -j 8 'awk '"'"'BEGIN { OFS="\t"; } {if($10=="CpG") print $1, $2, $2+1, $8;}'"'"' {} > {.}_CpG.bg' ::: rrbs_data/*.bed
LC_COLLATE=C

parallel -j 8 '[[ ! -f {.}_CpG.bg ]] && awk '"'"'BEGIN { OFS="\t"; } {if($10=="CpG") print $1, $2, $2+1, $8;}'"'"' {} | sort -k1,1 -k2,2n > {.}_CpG.bg' ::: rrbs_data/*.bed
parallel -j 8 '[[ ! -f {.}.bw ]] && bedGraphToBigWig {.}.bg hg19_chrom.sizes {.}.bw' ::: rrbs_data/*_CpG.bg


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
