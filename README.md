# Genetic implication of prenatal GABAergic and cholinergic neuron development in susceptibility to schizophrenia (2024)

This project was carried out in the Division of Psychological Medicine and Clinical Neurosciences (DPMCN). The paper is [here]() ADD WHEN PUBLISHED. 

***

## **snATACseq Data**

snATACseq data for this study are available in the following [EGA repository]():

Details on the snRNAseq data  / analyses for this study can be found [here](https://github.com/Dazcam/cameron_schizophr_bull_2023_snRNAseq).

***

## **GWAS Data**

See the following papers for GWAS data access:

+ [Schizophrenia](https://figshare.com/ndownloader/files/28169757)
+ [Autism](https://figshare.com/ndownloader/files/28169292)
+ [Major Depressive Disorder]() - Permission required at time of access
+ [ADHD](https://figshare.com/ndownloader/files/40036684)
+ [Bipolar Disorder]() - Permission required at time of access
+ [Height](https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz)

***

**Scripts**

1. [snRNAseq_GE_get_shi_data.sh](scripts/snRNAseq_GE_get_shi_data.sh) - Download Shi data
2. [snRNAseq_GE_prep_shi_data_for_Seurat.R](scripts/snRNAseq_GE_prep_shi_data_for_Seurat.R) - Clean and prepare Shi data
3. [snRNAseq_GE_shi_seurat.R](scripts/snRNAseq_GE_shi_seurat.R) - Generate Seurat objects (with batch correction) for cluster levels 1 and 2
4. [snRNAseq_GE_prep_enrich_test_files.R](scripts/snRNAseq_GE_prep_enrich_test_files.R) - Generate specificity scores and top 10% genesets for enrichment testing



***

Papers for public datasets used for confirmation of our cell assignments 

+ [Shi et al. 2021](https://www.science.org/doi/10.1126/science.abj6641) - snATACseq from fetal brain
 
***

Main scripts of interest

Alignment (Cell Ranger ATAC) and basic snATACseq processing (ArchR)

1. [snATACseq_cellRanger.smk](workflow/rules/snATACseq_cellRanger.smk) - Run Cell Ranger `atac-count` on fastQ files 
2. [snATACseq_archR_processing.smk](workflow/rules/snATACseq_archR_processing.smk) - Run ArchR QC and batch correction
    + [snATACseq_QC.R](workflow/scripts/snATACseq_QC.R) - Run ArchR QC 
3. [snATACseq_local_testing.R](workflow/scripts/snATACseq_local_testing.R) - Run main ArchR snATACseq processes locally

Overlap of fine-mapped SCZ-assoc. SNPs in co-accessible OCRs 

4. [snATACseq_map_PGC3_SCZ_finemapped_SNPs_to_peaks.R](workflow/scripts/snATACseq_map_PGC3_SCZ_finemapped_SNPs_to_peaks.R) - Map SCZ SNPs to cell specific OCRs
5. [snATACseq_annotate_archR_coaccesible_peaks.R](workflow/scripts/snATACseq_annotate_archR_coaccesible_peaks.R) - Annotate cA peaks and pull out peak pairs containing SCZ SNP
6. [snATACseq_find_overlapping_peaks.R](workflow/scripts/snATACseq_find_overlapping_peaks.R) - Compare OCRs to publicly available GE OCRs and run motif analysis

SLDSR

7. [snATACseq_LDSR.smk](workflow/rules/snATACseq_LDSR.smk) - Run sLDSR integrating cell specific OCRs with SCZ and Height GWAS data

