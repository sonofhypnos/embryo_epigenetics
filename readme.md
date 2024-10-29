# Next steps

- [ ] add documentation on data and source files
- [ ] At the end we probably want to export all of the data to an amazon bucket and then delete the SSD here. Unfortunately I did not document all of the steps of things to install (like all of the things I needed for bismark (I'll try to roughly list some of them here, so I can later do this in a script))
  - [ ] manually install igzip (./autoconf) (This worked, but there were like 4 steps inbetween where it would complain about some dependency not being installed, but installing them through apt just worked without any hiccups.)
  - [ ] install conda (install bioconda etc) (install whatever file there was to document conda installation)
  - [ ] install conda, make an new environment and run:
```
conda install -c bioconda bismark bedtools trim-galore sra-tools
```
  - [ ] install bismark (through conda, ask claude/look in your claude chat)
  - [ ] I also installed my dotfiles and other things for my ease of use (very optional step)
- [ ] validate that the methylation regions I used with the bed files make sense 
- [ ] I generated test files like this:
```
seqtk sample -s100 ../fastq/SRR5720828_1.fastq.gz 10000 > test.fastq   
(epi_env) ➜  methylation_analysis_1.0 git:(main) ✗ mv test.fastq ../fastq/test_1.fastq
(epi_env) ➜  methylation_analysis_1.0 git:(main) ✗ seqtk sample -s100 ../fastq/SRR5720828_2.fastq.gz 10000 > ../fastq/test_2.fastq
```


