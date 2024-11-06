for cell_type in Zygote "2cell" "4cell" "8cell" Morula ICM Sperm Oocyte TE hESC "alpha-Amanitin"; do
    for meth_type in "ACG.TCG" "GCA.GCC.GCT"; do
	    echo converting $cell_type and methylation: $meth_type...
	    ./wig_to_big_wig.sh "${cell_type}_${meth_type}.wig" hg19.chrom.sizes "${cell_type}_${meth_type}.bw"
	    echo done converting $cell_type and methylation: $meth_type
    done
done
