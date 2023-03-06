# Fetal brain snATACseq analyses (version 3)

Focuusing on GE snATACseq data (from V1) using ArchR. 

***

Papers for public datasets used for confirmation of our cell assignments 

+ [Shi et al. 2021](https://www.science.org/doi/10.1126/science.abj6641) - snATACseq from fetal brain
 
***

Main scripts of interest

ArchR

1. [snATACseq_cellRanger.smk](workflow/rules/snATACseq_cellRanger.smk) - Run Cell Ranger `atac-count` on fastQ files 

