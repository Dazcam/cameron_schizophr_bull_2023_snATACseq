#--------------------------------------------------------------------------------------
#
#     Map PGC3 SCZ fine mapped SNPs to snATACseq peaks
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# 1. Download PGC3 SCZ GWAS fine mapped SNPs
#    16 SNPs with 1:28690628_T_C encoding removed for now
# 2. Map SNPs to hg38 using BioMart
# 3. Check for overlap of PGC3 SCZ finemapped SNPs in snATACseq peaks (250 bp ext)
# 4. Create binary df for whether SNP is in/not in peak for all cell types
# 5. Run Wilcoxon rank sum test to test for difference in posterior probability means
#    between PGC3 SCZ SNPs in peaks v.s. rest of PGC3 SNPs not in peaks 

# Need to re-annotate and add the 23 SNPs with wrong encoding in PGC3

##  Load Packages  --------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(Repitools)
library(biomaRt)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
IN_DIR <- "~/Desktop/fetal_brain_snATACseq_V3_010323/resources/public_datasets/trubestskoy_2022/2020-08-14908C-s11/"
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_V3_010323/results/05PEAKS/"
SNP_DIR <- paste0(PEAK_DIR, 'finemapped_SNPs/')
PEAK_EXTENSION <- c('ext250bp')
CELL_TYPES <- c("CGE", "LGE", "MGE", "Progenitor")

##  Create directories  ----------------------------------------------------------------
dir.create(SNP_DIR)

##  Read in PGC3 index SNPs one is not an rsID just removed it  ------------------------
snps <- read_excel(paste0(IN_DIR, 'Supplementary Table 11.xlsx'), 
                   sheet = 'ST11a 95% Credible Sets') %>%
  dplyr::select(rsid) %>% # 20750 SNPs
  filter(!grepl(':|_', rsid)) %>% # 23 SNPs with 1:28690628_T_C encoding removed for now
  base::as.data.frame(snps) %>%
  distinct(rsid, .keep_all = TRUE) %>%
  arrange(rsid) %>%
  pull()

cat(paste0(length(snps), ' SNPs retained from PGC3 table.\n'))

snps_norsID <- read_excel(paste0(IN_DIR, 'Supplementary Table 11.xlsx'),
                          sheet = 'ST11a 95% Credible Sets') %>%
   dplyr::select(rsid) %>%
   filter(grepl(':|_', rsid))

 cat(paste0(length(snps_norsID %>% pull), ' SNPs with no rsIDs.\n'))

# # Get rsIDs for snps_norsID
# ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# 
# ## Create an example set of coordinates as a data frame
# snps_coords <- snps_norsID %>%
#   separate(rsid, c('Chr', 'Base'), sep = ':') %>%
#   separate(Base, c('Start', 'Alleles'), sep = '_', extra = 'merge') %>%
#   dplyr::select(-Alleles) %>%
#   mutate(Stop = Start) %>%
#   apply(1, paste, collapse = ":")
# 
# # Very slow some SNPs not found
# snps_extra[,5] <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
#                     filters = 'chromosomal_region', 
#                     values = snps_coords[5], 
#                     mart = ensembl)  

## Get hg38 base postions for rsIDs using biomaRt -------------------------------------
## ~15-40 mins per 50K SNPs - note sometimes crashes when BiomaRt in heavy use
cat('\nUsing BiomaRt to get hg38 base postions for SNP rsIDs  ... \n')
ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# Need to batch the query: https://support.bioconductor.org/p/23684/
SNPs <- getBM(attributes=c("refsnp_id",
                           "chr_name",
                           "chrom_start",
                           "chrom_end"),
                filters = "snp_filter", 
                values = snps, 
                mart = ensembl, 
                uniqueRows = TRUE)

cat(paste0(nrow(SNPs), ' SNPs retained after lift over to hg38.\n'))

# Some SNPs include duplicates due to chr patches - I removed these
cat('\nRemoving SNPs on CHR patches ... \n')
snps_no_patches <- SNPs %>%
  filter(!grepl('_', chr_name)) %>%
  dplyr::select(-chrom_end) %>%
  dplyr::rename('rsid' = refsnp_id, 'hg38_base_position' = chrom_start)
cat(paste0(nrow(snps_no_patches), ' SNPs retained. \n'))

# Add posterior probability and index SNPs columns 
SNPs_join <- snps_no_patches %>% 
  left_join(read_excel(paste0(IN_DIR, 'Supplementary Table 11.xlsx'), 
                       sheet = 'ST11a 95% Credible Sets')) %>%
  dplyr::select(rsid, chr_name, hg38_base_position, index_snp, 
                finemap_posterior_probability)

# Write to file
write_tsv(as.data.frame(SNPs_join), paste0(SNP_DIR, 'pgc3_scz_finemapped_SNPs_hg38.tsv'))


## Check for overlap of SNPs in snATACseq peaks of individual cell types  -------------
# Could speed be improved using Granges? ChipSeeker? findOverlaps? ccpR?
for (CELL_TYPE in CELL_TYPES) {
  
    for (EXT in PEAK_EXTENSION) {
    
      cat(paste0('\nLoading peaks for ', CELL_TYPE, ' ... \n'))
      peaks_df <- read_tsv(paste0(PEAK_DIR, CELL_TYPE,'.hg38.', EXT, '.bed'), col_names = FALSE)
      colnames(peaks_df) <- c('chr', 'start', 'end', 'name', 'score', 'strand')
    
      cell_overlaps <- data.frame()
      
      cat(paste0('\nChecking for SNP overlaps in ', CELL_TYPE, ' ... \n'))
      for (i in 1:nrow(snps_no_patches)) {
        
        BASE_POSITION <- snps_no_patches$hg38_base_position[i]
        CHR <- snps_no_patches$chr_name[i]
        cat(paste0('\nSNP: ', 
                   snps_no_patches$rsid[i], ', position: ',
                   BASE_POSITION, ' ... \n'))
        
        overlaps <- filter(peaks_df, start <= BASE_POSITION, end >= BASE_POSITION, chr == paste0('chr', CHR))
        
        if (nrow(overlaps) > 0) { 
          
          overlaps <- cbind(overlaps, snps_no_patches[i,])
          print(overlaps) 
          
        }
        
        cell_overlaps <- rbind(cell_overlaps, overlaps)
        
      } 
      
      cat(paste0('\nAll ', nrow(snps_no_patches), ' SNPs checked in ', CELL_TYPE, ' with peak ', EXT, '... \n'))
      
      # Add posterior probability and index SNP info to output
      cell_overlaps <- cell_overlaps %>% 
        left_join(snps <- read_excel(paste0(IN_DIR, 'Supplementary Table 11.xlsx'), 
                                     sheet = 'ST11a 95% Credible Sets')) %>%
        dplyr::select(rsid, chr, start, end, hg38_base_position, index_snp, finemap_posterior_probability)
      
      cat(paste0('\nWriting overlapping SNPs to file ... \n'))
      write_tsv(cell_overlaps, paste0(SNP_DIR, CELL_TYPE, 
                                      '_PGC3_SCZ_finemapped_SNP_peak_overlaps_', EXT , '.tsv'))
      assign(paste0(CELL_TYPE, '_', EXT, '_SNP_overlaps'), cell_overlaps)
    
    }
  
}


## Create binary count table to have all SNP/peak overlaps in one df  -----------------
# Create key of all SNPs overlapping peaks
for (EXT in PEAK_EXTENSION) {
  
  if (EXT == 'ext250bp') {
    
    cat('\nCreating key dataframe for all overlapping SNPs/peaks', EXT ,'... \n')
    all_SNPs_key_ext250bp_df <- rbind(CGE_ext250bp_SNP_overlaps, LGE_ext250bp_SNP_overlaps, 
                                 MGE_ext250bp_SNP_overlaps, Progenitor_ext250bp_SNP_overlaps)
    all_SNPs_key_ext250bp <- all_SNPs_key_ext250bp_df %>%
      arrange(rsid) %>%
      distinct(rsid)
  
  } else {
    
    cat('\nCreating key dataframe for all overlapping SNPs/peaks', EXT ,'... \n')
    all_SNPs_key_ext250bp_df <- rbind(CGE_ext500bp_SNP_overlaps, LGE_ext500bp_SNP_overlaps, 
                                 MGE_ext500bp_SNP_overlaps, Progenitor_ext500bp_SNP_overlaps)
    all_SNPs_key_ext250bp <- all_SNPs_key_500_df %>%
      arrange(rsid) %>%
      distinct(rsid)
    
  }
  
}


# Create binary df - SNP in/not in peak for each cell
cat('\nCreating binary df for whether SNP is in/not in peak ... \n')
for (EXT in c(PEAK_EXTENSION)) {
  
  # Get key def for extension
  SNPS_KEY_DF <- get(paste0('all_SNPs_key_', EXT))
  
  # Needed to reset df for each extension
  if (exists('all_SNPs_binary_df')) { rm(all_SNPs_binary_df) }
  
  for (CELL_TYPE in CELL_TYPES) {
  
  cat('\nObtaining binary counts for:', CELL_TYPE, 'peak', EXT, '... \n')
  # Test vector of rsIDs for cell type
  snp_test <- get(paste0(CELL_TYPE, '_', EXT, '_SNP_overlaps')) %>%
    pull(rsid)
  
    if (exists('all_SNPs_binary_df')) {
      
      all_SNPs_binary_df <- all_SNPs_binary_df %>%
        rowwise() %>%
        mutate(!!CELL_TYPE := ifelse(rsid %in% snp_test, 1, 0)) %>%
        ungroup()
    
    } else {
      
      all_SNPs_binary_df <- SNPS_KEY_DF %>%
        rowwise() %>%
        mutate(!!CELL_TYPE := ifelse(rsid %in% snp_test, 1, 0)) %>%
        ungroup()
      
    }
  
  }
  
  cat('\nAssigning binary counts df for peak ext', EXT, '... \n')
  assign(paste0('all_SNPs_binary_', EXT, '_df'), all_SNPs_binary_df)
  
}

# Add index SNP and posterior probability (PP) data to binary dfs and run 
# Wilcoxon tests comparing difference in PPs of SNPs in peaks compared
# to SNPs not in peaks
for (EXT in c(PEAK_EXTENSION)) {
  
  SNPS_BINARY_DF <- get(paste0('all_SNPs_binary_', EXT, '_df'))
  SNPS_KEY_DF <- get(paste0('all_SNPs_key_', EXT, '_df'))
  
  # Add index SNP and posterior probability (PP) data to binary dfs
  SNPS_FINAL_DF_NO_PEAKS <- SNPS_BINARY_DF %>%
    left_join(SNPS_KEY_DF) %>%
    dplyr::select(-start, -end) %>%
    relocate(rsid, hg38_base_position, index_snp, finemap_posterior_probability) %>%
    distinct()
  
  SNPS_FINAL_DF_WITH_PEAKS <- SNPS_BINARY_DF %>%
    left_join(SNPS_KEY_DF) %>%
    relocate(rsid, hg38_base_position, index_snp, finemap_posterior_probability) %>%
    distinct()
  
  # Write binary dfs
  cat('\nWriting binary count tables ... \n')
  write_tsv(SNPS_FINAL_DF_NO_PEAKS, 
            paste0(SNP_DIR, 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_', 
                   EXT, '_SNPs_only.tsv'))
  
  write_tsv(SNPS_FINAL_DF_WITH_PEAKS, 
            paste0(SNP_DIR, 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_', 
                   EXT, '_with_peaks.tsv'))
  
  ## Run 2-sided Wilcoxon tests
  SNPS_FINAL_DF_rsIDs <- SNPS_FINAL_DF_NO_PEAKS %>%
    pull(rsid)
  
  SNPS_FINAL_DF_PPs <- SNPS_FINAL_DF_NO_PEAKS %>%
    pull(finemap_posterior_probability)
  
  PGC3_SNPS_PPs <- read_excel(paste0(IN_DIR, 'Supplementary Table 11.xlsx'), 
                              sheet = 'ST11a 95% Credible Sets') %>%
    dplyr::select(rsid, finemap_posterior_probability) %>%
    filter(!grepl(':|_', rsid)) %>% # 16 SNPs with 1:28690628_T_C encoding removed for now
    base::as.data.frame(snps) %>%
    distinct(rsid, .keep_all = TRUE) %>%
    arrange(rsid) %>%
    filter(!rsid %in% SNPS_FINAL_DF_rsIDs) %>%
    pull(finemap_posterior_probability)
  
  cat('\n\nWilcox test for SNPs ext:', EXT, '\n\n')
  
  print(wilcox.test(SNPS_FINAL_DF_PPs, PGC3_SNPS_PPs, alternative = 't'))
  
  cat('\n\n')
  
}

cat('Done.')


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------