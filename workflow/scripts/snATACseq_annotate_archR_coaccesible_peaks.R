#--------------------------------------------------------------------------------------
#
#     Annotate ArchR correlated co-accesible peaks
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# 1. Read in snATACseq cA peaks and metadata
# 2. Combine all peak information for cA correlated peaks in single df
# 3. Find cA peak pairs where one peak contains PGC3 SCZ finemapped SNP
# 4. Annotate cA peak pairs containing SNPs
# 5. Extract correlated peaks where peak NOT containing SNP are annotated as promoter

# _cA_peaks_df - queryHits and subjectHits cols denote index of the two correlated peaks
# _cA_peaks_metadata - indexes of queryHits and subjectHits mentioned above apply to this


##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86) #hg38
library(org.Hs.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(Repitools) # Granges to df
#library(diffloop) # Add/rm chr to granges seqlevels

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_V3_010323/results/05PEAKS/"
SNP_DIR <- paste0(PEAK_DIR, "finemapped_SNPs/")
REGIONS <- "GE"

##  Load data  ----------------------------------------------------------------
cat('\nLoading peak and SNP data ... \n')
for (REGION in REGIONS) {
  
  PEAKS <- readRDS(paste0(PEAK_DIR, REGION, '_cA_peaks_df.rds'))
  METADATA <- readRDS(paste0(PEAK_DIR, REGION, '_cA_peaks_metadata.rds'))
  
  assign(paste0(REGION, '_cA_peaks'), PEAKS)
  assign(paste0(REGION, '_cA_metadata'), METADATA)
    
}

SNPS <- read_tsv(paste0(SNP_DIR, 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext250bp_SNPs_only.tsv'))
SNPS_WITH_PEAKS <- read_tsv(paste0(SNP_DIR, 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext250bp_with_peaks.tsv'))

##  Combine info on cA peaks and metadata dfs  ----------------------------------------
# Get peak indices for query and subject hits in index files - these may be the same - check
cat('\n\nCombine peak info stored in cA peaks and metadata dfs ... \n')
for (REGION in REGIONS) {
  
  cat('\nRunning ', REGION, ' ...\n')
  
  CA_PEAKS <- get(paste0(REGION, '_cA_peaks'))
  CA_METADATA <- get(paste0(REGION, '_cA_metadata'))

  # Get indices for query and subject peaks
  query_peaks_idx <- CA_PEAKS$queryHits
  subject_peaks_idx <- CA_PEAKS$subjectHits
  
  # Pull peak information for indices from metadata
  query_peaks_info <- CA_METADATA[query_peaks_idx] 
  subject_peaks_info <- CA_METADATA[subject_peaks_idx] 
  
  # Remove rownames
  names(query_peaks_info) <- NULL
  names(subject_peaks_info) <- NULL
  
  # Convert peak info to df
  query_peaks_info_no_chr <- annoGR2DF(query_peaks_info) %>%
    dplyr::select(start, end)
  subject_peaks_info_no_chr <- annoGR2DF(subject_peaks_info) %>%
    dplyr::select(start, end)
  colnames(query_peaks_info_no_chr) <- c("query_start", "query_end")
  colnames(subject_peaks_info_no_chr) <- c("subject_start", "subject_end")
  
  # Save peak info with chr
  query_peaks_info_with_chr <- annoGR2DF(query_peaks_info) %>%
    dplyr::select(chr, start, end)
  subject_peaks_info_with_chr <- annoGR2DF(subject_peaks_info) %>%
    dplyr::select(chr, start, end)
  
  # Convert peaks S4 df to S3
  CA_PEAKS <- as.data.frame(CA_PEAKS)
  CA_PEAKS_NO_CORR <- cbind(query_peaks_info_with_chr, subject_peaks_info_with_chr)
  
  # Check if chr cols identical
  identical(CA_PEAKS_NO_CORR[,1], CA_PEAKS_NO_CORR[,4])
  colnames(CA_PEAKS_NO_CORR ) <- c("chr", "query_start", "query_end", "chr_ext", 
                          "subject_start", "subject_end")
  CA_PEAKS_NO_CORR  <- CA_PEAKS_NO_CORR  %>%
    select(-chr_ext)
  
  CA_PEAKS_WITH_CORR <- cbind(query_peaks_info_no_chr, subject_peaks_info_no_chr, CA_PEAKS) 
  
  cat('Total cA peak pairs in', REGION, 'after conversion:', nrow(CA_PEAKS), '\n')
  
  assign(paste0(REGION, '_cA_peaks_combined'), as_tibble(CA_PEAKS_NO_CORR))
  assign(paste0(REGION, '_cA_peaks_combined_with_cor'), as_tibble(CA_PEAKS_WITH_CORR))
  assign(paste0(REGION, '_cA_query_peaks_combined'), as_tibble(query_peaks_info_with_chr))
  assign(paste0(REGION, '_cA_subject_peaks_combined'), as_tibble(subject_peaks_info_with_chr))
  
  write_tsv(CA_PEAKS, paste0(PEAK_DIR, REGION, '_cA_peaks_combined.tsv'))
  
  # Clean up
  rm(list = ls(pattern="^CA_"))
  rm(METADATA, PEAKS)
  rm(list = ls(pattern="^subject_"))
  rm(list = ls(pattern="^query_"))
  
}

##  Find cA peak pairs where one peak contains PGC3 SCZ finemapped SNP  ---------------
cat('\n\nFind cA peak pairs where one peak contains PGC3 SCZ finemapped SNP ... \n')
for (REGION in REGIONS) {
  
  # Note there may be a fair bit of redundancy here
  # sum(sort(subject_peaks_idx) == sort(query_peaks_idx)) = 153886
  # length(query_peaks_idx) = 153894
  
  # Pull out ranges from the ArchR cA metadata dataframe
  ALL_CA_PEAKS <- get(paste0(REGION, '_cA_metadata'))
  PEAKS_COMBINED <- get(paste0(REGION, '_cA_peaks_combined_with_cor')) %>%
    select(-Variability1, -Variability2, -TStat, Pval, -VarQuantile1, -VarQuantile2) %>%
    rename(chr = seqnames)
  names(ALL_CA_PEAKS) <- NULL
  
  # Create GRanges object for all PGC3 SCZ SNPs that overlap snATACseq peak
  # use only hg38 base position - checked if all BPs are unique in SNPs
  cat('\n\nCreate GRanges object for finemapped SNPs ... \n')
  SNPS_add_range <- SNPS %>%
    mutate(end = hg38_base_position) %>%
    rename(start = hg38_base_position) %>%
    relocate(rsid, chr, start, end)
  SNPS_granges <- makeGRangesFromDataFrame(SNPS_add_range, keep.extra.columns = FALSE,
                                           ignore.strand = TRUE)
  
  # Find overlaps using IRanges
  OVERLAPS <- findOverlaps(ALL_CA_PEAKS, SNPS_granges)
  
  # Combine overlapping peak and SNP info in single df - only base position at this stage 
  METADATA_OVERLAPS <- as.data.frame(ALL_CA_PEAKS[queryHits(OVERLAPS)]) %>%
    select(seqnames, start, end)
  SNP_OVERLAPS <- as.data.frame(SNPS_granges[subjectHits(OVERLAPS)]) %>%
    select(seqnames, start) %>%
    rename(hg38_base_position = start)
  OVERLAPS_DF <- cbind(METADATA_OVERLAPS, SNP_OVERLAPS) 
  colnames(OVERLAPS_DF) <- c("chr", "start", "end", "seqnames", "hg38_base_position")
  
  # Do chr columns match between SNPs and peaks?
  identical(as.vector(OVERLAPS_DF$chr), as.vector(OVERLAPS_DF$seqnames))

  # Rm duplicate chr column
  OVERLAPS_DF <- OVERLAPS_DF %>%
    select(-seqnames)
  
  # Double check overlaps
  for (SNP in 1:nrow(OVERLAPS_DF)) {
    
    test <- ifelse(OVERLAPS_DF$hg38_base_position[SNP] >= OVERLAPS_DF$start[SNP] & OVERLAPS_DF$hg38_base_position[SNP] <= OVERLAPS_DF$end[SNP], TRUE, FALSE) 
    print(test)
    
  }
  
  # Join with rsIDs - increase in entries here as multiple SNPs in single peaks
  OVERLAPS_DF_JOIN <- OVERLAPS_DF %>% left_join(SNPS) %>%
    select(chr, start, end, hg38_base_position, rsid)
  
  # Unique SNPs in final df
  cat('\nNumber of unique rsIDs after overlap of SNPs with cA metadata:',
      length(unique(OVERLAPS_DF_JOIN$rsid)), '\n')
 
  # Pull out peaks from query and subject in combined metadata peaks
  QUERY_JOIN <- OVERLAPS_DF_JOIN %>% 
    inner_join(PEAKS_COMBINED, by = c('chr' = 'chr', 'start' = 'query_start', 'end' = 'query_end'))
  
  SUBJECT_JOIN <- OVERLAPS_DF_JOIN %>% 
    inner_join(PEAKS_COMBINED, by = c('chr' = 'chr', 'start' = 'subject_start', 'end' = 'subject_end'))
  
  # Note that the QUERY_JOIN and SUBJECT_JOIN objects are the same - the metadata object must contain
  # the peak information both ways (this was before I added correlation values)
  identical(QUERY_JOIN[,1:7], SUBJECT_JOIN[,1:7] %>% rename(subject_start = query_start, subject_end = query_end))

  # Unique SNPs after joining with cA peaks
  cat('\nNumber of unique rsIDs after overlap of SNPs with cA metadata query:',
      length(unique(QUERY_JOIN$rsid)), '\n')
  
  assign(paste0(REGION, '_cA_peak_pairs_with_SNP'), as_tibble(QUERY_JOIN))

}

##  Annotate cA peak pairs containing SNPs  -------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cat('\n\nRunning regulatory peak annotation for cA peaks pairs containing SNPs ... \n')
for (REGION in REGIONS) {
  
  cat('\nRunning ', REGION, ' ...\n')
  
  # Load peaks - need to remove chr to match edb encoding for annotatePeaks
  cA_PEAKS <- get(paste0(REGION, '_cA_peak_pairs_with_SNP')) #%>%
   # mutate(chr = gsub("chr", "", chr))
  
  # Split df to run each set of cA peaks separately - peak with SNP is query set
  cat('\n\nSplit df to run each set of cA peaks separately ... \n')
  cA_PEAKS_with_SNP <- cA_PEAKS %>% 
    select(-subject_start, -subject_end)
  cA_PEAKS_no_SNP <- cA_PEAKS %>% 
    select(-start, -end) %>%
    rename(start = subject_start,
           end = subject_end)
    
  # Create GRanges object for each peak set
  cat('\n\nCreate GRanges object for each peak set ... \n')
  cA_PEAKS_with_SNP_GR <- makeGRangesFromDataFrame(cA_PEAKS_with_SNP, keep.extra.columns = TRUE,
                                                   ignore.strand = TRUE)
  cA_PEAKS_no_SNP_GR <- makeGRangesFromDataFrame(cA_PEAKS_no_SNP, keep.extra.columns = TRUE,
                                                  ignore.strand = TRUE)
  
  # Annotate peaks
  cat('\n\nAnnotating peaks with SNP ... \n\n\n')
  cA_PEAKS_with_SNP_ANN <- annotatePeak(cA_PEAKS_with_SNP_GR, tssRegion=c(-1000, 100),
                                  TxDb = txdb, annoDb = "org.Hs.eg.db", level = 'gene')
  print(cA_PEAKS_with_SNP_ANN)
  
  cat('\n\nAnnotating peaks without SNP ... \n\n')
  cA_PEAKS_no_SNP_ANN <- annotatePeak(cA_PEAKS_no_SNP_GR, tssRegion=c(-1000, 100),
                                  TxDb = txdb, annoDb = "org.Hs.eg.db", level = 'gene')

  cA_PEAKS_with_SNP_ANN_DF  <- as.data.frame(cA_PEAKS_with_SNP_ANN@anno) %>%
    mutate(rsid = cA_PEAKS$rsid,
           hg38_base_position = cA_PEAKS$hg38_base_position) %>%
    relocate(seqnames, start, end, hg38_base_position, rsid) %>%
    select(-strand, -width, -queryHits, -subjectHits) 
  
  cA_PEAKS_no_SNP_ANN_DF  <- as.data.frame(cA_PEAKS_no_SNP_ANN@anno) %>%
    mutate(rsid = cA_PEAKS$rsid,
           hg38_base_position = cA_PEAKS$hg38_base_position) %>%
    relocate(seqnames, start, end, hg38_base_position, rsid) %>%
    select(-strand, -width, -queryHits, -subjectHits)
  
  # Add SNPs info back into df
  assign(paste0(REGION, '_cA_peak_pairs_with_SNP_ann'), cA_PEAKS_with_SNP_ANN_DF)
  assign(paste0(REGION, '_cA_peak_pairs_no_SNP_ann'), cA_PEAKS_no_SNP_ANN_DF)
  
  #write_tsv(cA_PEAKS_with_SNP_ANN_DF, paste0(SNP_DIR, REGION, '_cA_regulatory_anns_with_finemapped_SNP.tsv'))
  #write_tsv(cA_PEAKS_no_SNP_ANN_DF, paste0(PEAK_DIR, REGION, '_cA_regulatory_anns_no_finemapped_SNP.tsv'))
  
}

# Generate final table
final_cA_table <- GE_cA_peak_pairs_with_SNP_ann %>%
  as_tibble() %>%
  select(-geneChr, -geneStart, -geneEnd, -geneLength, -geneId, -hg38_base_position,
         -distanceToTSS, -GENENAME, -geneStrand, -SYMBOL, -ENSEMBL) %>%
  rename(ann_pk_w_SNP = annotation,
         start_snp = start,
         end_snp = end) %>%
  left_join(GE_cA_peak_pairs_no_SNP_ann) %>%
  select(-geneChr, -geneStart, -geneEnd, -geneLength, -geneId, -hg38_base_position,
         -distanceToTSS, -GENENAME, -geneStrand, -Pval) %>%
  rename(ann_pk_no_SNP = annotation,
         start_no_snp = start,
         end_no_snp = end,
         chr = seqnames) %>%
  relocate(chr, start_snp, end_snp, start_no_snp, end_no_snp, rsid) %>%
  filter(ann_pk_no_SNP == 'Promoter') %>%
  arrange(FDR) 
  
write_tsv(final_cA_table, paste0(SNP_DIR, REGION, '_cor_cA_peaks_with_finemapped_SNP.tsv'))


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
