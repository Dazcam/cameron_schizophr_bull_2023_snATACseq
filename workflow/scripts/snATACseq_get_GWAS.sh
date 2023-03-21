#!/bin/bash

# SCZ data wave 3 v2 file is here: https://doi.org/10.6084/m9.figshare.14672178 - N taken as 161405 = European (EUR) and East Asian (ASN) ancestry cohorts only, see frq cols in sumstats
# Height data Wood et al + UK BioBank is here: https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files

# Set variables
SCZ_LOC=https://figshare.com/ndownloader/files/28169757
HEIGHT_LOC=https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz
OUTDIR=../resources/gwas/

mkdir ${OUT_DIR}

# Download GWAS
wget ${SCZ_LOC} -P ${OUTDIR}
wget ${HEIGHT_LOC} -P ${OUTDIR}  

# Unpack GWAS
mv ${OUTDIR}28169757 ${OUTDIR}PGC3_SCZ_wave3_public.v2.tsv.gz
zcat ${OUTDIR}Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz > ${OUTDIR}Meta-analysis_Wood_et_al+UKBiobank_2018.txt
zcat ${OUTDIR}PGC3_SCZ_wave3_public.v2.tsv.gz > ${OUTDIR}PGC3_SCZ_wave3_public.v2.tsv
 
 
