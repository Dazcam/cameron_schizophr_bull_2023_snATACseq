#!/bin/bash

# Run from workflow directory

# Set variables
LDSR_ROOT=../resources/ldsr/
LDSR_DIR=../resources/ldsr/reference_files

# Clone repo
git clone git@github.com:bulik/ldsc.git ${LDSR_ROOT}

# Create dir for LDSR support files
mkdir -p ${LDSR_DIR}

# Download support files
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baseline_v1.2_ldscores.tgz -P ${LDSR_DIR}
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_plinkfiles.tgz -P ${LDSR_DIR}
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/w_hm3.snplist.bz2 -P ${LDSR_DIR}
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/weights_hm3_no_hla.tgz -P ${LDSR_DIR}
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_frq.tgz -P ${LDSR_DIR}
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/hapmap3_snps.tgz -P ${LDSR_DIR}
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 -P ${LDSR_DIR}

# Unpack support files
tar zxf ${LDSR_DIR}/1000G_Phase3_baseline_v1.2_ldscores.tgz -C ${LDSR_DIR}/
tar zxf	${LDSR_DIR}/1000G_Phase3_plinkfiles.tgz -C ${LDSR_DIR}
tar zxf	${LDSR_DIR}/weights_hm3_no_hla.tgz -C ${LDSR_DIR}
tar zxf	${LDSR_DIR}/1000G_Phase3_frq.tgz -C ${LDSR_DIR}
tar zxf ${LDSR_DIR}/hapmap3_snps.tgz -C ${LDSR_DIR}
tar xf  ${LDSR_DIR}/eur_w_ld_chr.tar.bz2 -C ${LDSR_DIR}
bzip2 -dkc ${LDSR_DIR}/w_hm3.snplist.bz2 > ${LDSR_DIR}/w_hm3.snplist

# Rename baseline files dir
mv ${LDSR_DIR}/baseline_v1.2/ ${LDSR_DIR}/baseline_v1.2_1000G_Phase3

# rm originals 
rm -rf ${LDSR_DIR}/1000G_Phase3_baseline_v1.2_ldscores.tgz
rm -rf ${LDSR_DIR}/1000G_Phase3_plinkfiles.tgz
rm -rf ${LDSR_DIR}/weights_hm3_no_hla.tgz
rm -rf ${LDSR_DIR}/1000G_Phase3_frq.tgz
rm -rf ${LDSR_DIR}/hapmap3_snps.tgz
rm -rf ${LDSR_DIR}/eur_w_ld_chr.tar.bz2
rm -rf ${LDSR_DIR}/w_hm3.snplist.bz2

# Need to modify snplist - see here https://github.com/bulik/ldsc/issues/96
cut -f1 ${LDSR_DIR}/w_hm3.snplist > ${LDSR_DIR}/w_hm3.snplist_rsIds
