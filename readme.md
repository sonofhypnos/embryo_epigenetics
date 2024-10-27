# Next steps

- [ ] At the end we probably want to export all of the data to an amazon bucket and then delete the SSD here. Unfortunately I did not document all of the steps of things to install (like all of the things I needed for bismark (I'll try to roughly list some of them here, so I can later do this in a script))
  - [ ] manually install igzip (./autoconf) (This worked, but there were like 4 steps inbetween where it would complain about some dependency not being installed, but installing them through apt just worked without any hiccups.)
  - [ ] install conda (install bioconda etc)
  - [ ] install bismark (through conda, ask claude/look in your claude chat)
  - [ ] I also installed my dotfiles and other things for my ease of use (very optional step)
- [ ] validate that the methylation regions I used with the bed files make sense

