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

# Experiment I ran to compare correlation between different samples in the RRBS dataset

