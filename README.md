# Fetal brain snATACseq analyses (version 3) - ongoing

This project was carried out in the Division of Psychological Medicine and Clinical Neurosciences (DPMCN). (https://www.biologicalpsychiatryjournal.com/article/S0006-3223(22)01404-4/fulltext). The workflow follows the the snakemake [distribution and reproducibility](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) recommendations. 

***

A snakemake (and local) pipeline to process snATACseq data. Utilising the following packages:

+ [Snakemake 6.6.1](https://snakemake.readthedocs.io/en/stable/)
+ [CellRanger atac 2.1.0](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/algorithms/overview)
+ [ArchR 1.0.2](https://www.archrproject.com) 
+ [Signac 1.9.0](https://stuartlab.org/signac/)
+ [ArchR2Signac](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
+ [MAGMA 1.08](https://ctg.cncr.nl/software/magma)
+ [LD Score Regression 1.0.1](https://github.com/bulik/ldsc)
+ etc.

***

**Data**

Focuusing on GE snATACseq data (from V1) using ArchR. 

***

Papers for public datasets used for confirmation of our cell assignments 

+ [Shi et al. 2021](https://www.science.org/doi/10.1126/science.abj6641) - snATACseq from fetal brain
 
***

Main scripts of interest

ArchR

1. [snATACseq_cellRanger.smk](workflow/rules/snATACseq_cellRanger.smk) - Run Cell Ranger `atac-count` on fastQ files 

