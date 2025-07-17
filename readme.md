# Info


# Installation

Run
```
conda create -n epi_env python bismark bedtools trim-galore sra-tools wiggletools -c bioconda
```

  - [ ] Optional: manually install igzip (./autoconf) (This worked, but there were 4 steps inbetween where igzip complain about some dependency not being installed, but installing them through apt just worked without any hiccups.)

# Overview

- `src/wgbs` contains scripts to aggregate the files from the ['Single-cell multi-omics sequencing of human early embryos' dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100272) 
  Scripts under src/wgbs were originally written to be run in the same directory as the input files. If they are rerun from `src/wgbs` the paths need to be adjusted. The wgbs data should be saved under `/data`.
- `src/rrbs` contains the script to aggregate the rrbs data from [The DNA methylation landscape of human early embryos](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49828). The data should be saved under `/rrbs_data`
- `src/clustering.py` is the script to cluster the bisulfite sequencing data from both the whole genome and from the wgbs and the rrbs data as well as the rrbs data from the supersox paper that should be saved under `/data/supersox/rrbs/`.
- `src/cosine_similarity.py` runs cosine similarity between the wgbs, the rrbs and the supersox paper. New datasets can be added by adding an instance of the DatasetParams class to AllDatasetParams. It's then pretty easy to plot the dataset with the plotting functions (see the main function in `src/cosine_similarity.py`)

# Experiment I ran to compare correlation between different samples in the RRBS dataset

# Log
2025-03-19:
Using short shell scripts like this example to transform oocyte data to hg38 to compare it with Oocyte data from the imprintome paper (Note this script depends on some files like hg19 to be in a specific place):

```sh
# Convert WIG to bedGraph, sort, and perform liftOver
wiggletools write_bg - Oocyte_ACG.TCG.wig | \
sort -k1,1 -k2,2n > temp.bedGraph && \
liftOver temp.bedGraph ../hg19ToHg38.over.chain temp_lifted.bedGraph Oocyte_ACG.TCG_unlifted.txt && \
sort -k1,1 -k2,2n temp_lifted.bedGraph > temp_lifted_sorted.bedGraph

# Merge overlapping regions using bedtools
bedtools merge -i temp_lifted_sorted.bedGraph -c 4 -o mean > temp_lifted_merged.bedGraph

# Convert merged bedGraph to bigWig
bedGraphToBigWig temp_lifted_merged.bedGraph ../hg38.chrom.sizes Oocyte_ACG.TCG_hg38.bw

# Clean up temporary files
rm temp.bedGraph temp_lifted.bedGraph temp_lifted_sorted.bedGraph temp_lifted_merged.bedGraph
```

I also used off-set, so that 0 could be displayed properly for files, but then I noticed, that in igv, one can get the same behaviour by just setting the datarange to: min=-offset and mid=-offset and then the line that is usually at 0, will be at -offset.

# Other 

  I generated test files like this:
  ```
  seqtk sample -s100 ../fastq/SRR5720828_1.fastq.gz 10000 > test.fastq   
  (epi_env) ➜  methylation_analysis_1.0 git:(main) ✗ mv test.fastq ../fastq/test_1.fastq
  (epi_env) ➜  methylation_analysis_1.0 git:(main) ✗ seqtk sample -s100 ../fastq/SRR5720828_2.fastq.gz 10000 > ../fastq/test_2.fastq
  ```
