#!/bin/bash

# Process each cell type and methylation type separately
for cell_type in Sperm Oocyte Zygote "2cell" "4cell" "8cell" Morula ICM TE hESC "alpha-Amanitin"; do
    for meth_type in "ACG.TCG" "GCA.GCC.GCT"; do
        # Get all replicates for this condition
        files=($(ls *${cell_type}*${meth_type}*.bw))
        echo $files
        
        if [ ${#files[@]} -gt 0 ]; then
            echo "Processing ${cell_type} ${meth_type}"
	    echo "run command: wiggletools ratio sum ${files[@]} : sum map unit map offset 1 map default -1 ${files[@]} > ${cell_type}_${meth_type}.wig"
	    wiggletools ratio sum "${files[@]}" : sum map unit map offset 1 map default -1 "${files[@]}" > "${cell_type}_${meth_type}.wig"
	    wig2bed < "${cell_type}_${meth_type}.wig" > "${cell_type}_${meth_type}.bed"
        fi
    done
done

# Generate UCSC hub.txt
# (reuse our previous Python code for this part since it's just text generation)
