# Info

Scripts under src/wgbs were originally written to be run in the same directory as the input files. If they are rerun the paths need to be adjusted.

# Next steps

- [ ] add documentation on data and source files 
  - [ ] manually install igzip (./autoconf) (This worked, but there were like 4 steps inbetween where it would complain about some dependency not being installed, but installing them through apt just worked without any hiccups.)
  - [ ] install conda (install bioconda etc) (install whatever file there was to document conda installation)
  - [ ] install conda, make an new environment and run:
```
conda install -c bioconda bismark bedtools trim-galore sra-tools
```
  - [ ] install bismark (through conda, ask claude/look in your claude chat)
  - [ ] I also installed my dotfiles and other things for my ease of use (very optional step)
  - [ ] I generated test files like this:
  ```
  seqtk sample -s100 ../fastq/SRR5720828_1.fastq.gz 10000 > test.fastq   
  (epi_env) ➜  methylation_analysis_1.0 git:(main) ✗ mv test.fastq ../fastq/test_1.fastq
  (epi_env) ➜  methylation_analysis_1.0 git:(main) ✗ seqtk sample -s100 ../fastq/SRR5720828_2.fastq.gz 10000 > ../fastq/test_2.fastq
  ```

  - [ ] wiggletools was added as another dependency. wiggletools is really annoying to install later, because it won't find the correct c-libraries, but symlinking the wrong libraries to wiggletools seems to just work
- [ ] refactor data directories
- [ ] add documentation on where the most important processed data files are (like methylation data for oocytes)

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

