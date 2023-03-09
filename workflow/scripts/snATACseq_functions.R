
#library(docstring) # Required for functions to be searchable by ?

set_filter_params <- function(
    
  FILTER_CODE = NULL ) {
  
  #' Set hard filters for tss and min frags  
  #' 
  #' @description Set hard filters for tss and min frags. Options include:
  #'  hard_trim, ziffra_only_default, frag_3.5.
  #' 
  #' @param FILTER_CODE: Either hard_trim, ziffra_only_default, frag_3.5
  
  if (FILTER_CODE == 'hard_trim') {
    
    pass_510 <- which(archR$Sample == 'PFC_510' & archR$TSSEnrichment >= 10 & archR$log10nFrags >= 4)
    pass_611 <- which(archR$Sample == 'PFC_611' & archR$TSSEnrichment >= 7.5 & archR$log10nFrags >= 3.75)
    pass_993 <- which(archR$Sample == 'PFC_993' & archR$TSSEnrichment >= 8 & archR$log10nFrags >= 3.5)
    pass_G17 <- which(archR$Sample == 'CTX_G17' & archR$TSSEnrichment >= 12)
    pass_G18 <- which(archR$Sample == 'CTX_G18' & archR$TSSEnrichment >= 12 & archR$log10nFrags >= 4)
    pass_G20 <- which(archR$Sample == 'PFC_G20' & archR$TSSEnrichment >= 12 & archR$log10nFrags >= 3.5)
    pass_G21 <- which(archR$Sample == 'CTX_G21' & archR$TSSEnrichment >= 10 & archR$log10nFrags >= 3.25)
    pass_T31 <- which(archR$Sample == 'PFC_T31' & archR$TSSEnrichment >= 15)
    pass_T34 <- which(archR$Sample == 'PFC_T34' & archR$TSSEnrichment >= 10.5 & archR$log10nFrags >= 3.8)
    
    all_pass <- c(pass_510, pass_611, pass_993, pass_G17, pass_G18,
                  pass_G20, pass_G21, pass_T31, pass_T34)
    
    return(all_pass)
    
  } else if (FILTER_CODE == 'ziffra_only_default') {
    
    pass_G17 <- which(archR$Sample == 'CTX_G17')
    pass_G18 <- which(archR$Sample == 'CTX_G18')
    pass_G20 <- which(archR$Sample == 'PFC_G20')
    pass_G21 <- which(archR$Sample == 'CTX_G21')
    pass_T31 <- which(archR$Sample == 'PFC_T31')
    pass_T34 <- which(archR$Sample == 'PFC_T34')
    
    all_pass <- c(pass_G17, pass_G18, pass_G20, pass_G21, pass_T31, pass_T34)
    
    return(all_pass)
    
  } else if (FILTER_CODE == 'frag_3.5') {
    
    all_pass <- which(archR$log10nFrags >= 3.5)
    
    return(all_pass)
    
  } else {
    
    print('No filter run. Filter options include: hard_trim, ziffra_only_default, frag_3.5.')
    
  }

}

create_umap_batch_plots <- function(
    
    ARCHR_ID = NULL, 
    CLUSTERS_ID =  NULL,
    UMAP_ID = NULL ){
  
  #' Creates a group UMAP plot  
  #' 
  #' @description Creates a group UMAP plot for each batch variable
  #' in ArchR metadata. Requires variables in BATCH_VAR vector. 
  #' 
  #' @param ARCHR_ID: ArchR object
  #' @param CLUSTERS_ID: Name of Cluster object stored in ArchR object
  #' @param UMAP_ID: Name of UMAP object stored in ArchR object


  UMAP_batch_plot_list <- list()
  
  for (NUM in 1:length(BATCH_VARS)) {
    
    cat('\nPlotting UMAP for:', BATCH_VARS[NUM], '\n\n')
  
    PLOT <- plotEmbedding(
      ArchRProj = ARCHR_ID,
      colorBy = "cellColData",
      name = BATCH_VARS[NUM],
      embedding = UMAP_ID,
      plotAs = 'points',  # See streaking issue https://github.com/GreenleafLab/ArchR/issues/1731
    ) + ggtitle(BATCH_VARS[NUM])
    
    UMAP_batch_plot_list[[NUM]] <- PLOT
    
  }

  UMAP_batch_plot_list[[6]] <- plotEmbedding(ArchRProj = ARCHR_ID, colorBy = "cellColData", 
                                             name = CLUSTERS_ID, 
                                             embedding = UMAP_ID,
                                             plotAs = 'points') +
    Seurat::NoLegend() + ggtitle('Clusters')
  
  UMAP_batch_plot_list <- lapply(UMAP_batch_plot_list, function(x){
    
    x + theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        legend.title = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
      )
  })

  assign(paste0('UMAP_batch_plot_list_', CLUSTERS_ID), UMAP_batch_plot_list, .GlobalEnv)
  cat('\nReturned:', paste0('UMAP_batch_plot_list_', CLUSTERS_ID), '\n\n')
  
}

create_archR_group_plot <- function(
    
  ARCHR_ID = NULL, 
  CLUSTERS_ID =  NULL,
  UMAP_ID = NULL ){
  
  #' Creates a group QC plot  
  #' 
  #' @description Creates a QC group plot for archR containing UMAPs bu cluster 
  #' and sample, a stacked bar plot and confusion matrix showing the 
  #' contribution from each sample to each cluster 
  #' 
  #' @param ARCHR_ID: ArchR object
  #' @param CLUSTERS_ID: Name of Cluster object stored in ArchR object
  #' @param UMAP_ID: Name of UMAP object stored in ArchR object
  
  ## Clusters - reporting  -----------------------------------------------------
  # Cluster counts - after Iterative LSI based clustering
  
  cluster_vector <- unname(unlist(getCellColData(ARCHR_ID)[CLUSTERS_ID]))
  sample_vector <- unname(unlist(getCellColData(ARCHR_ID)['Sample']))
  
  cat('\nCreating tables and plots for Iterative LSI based clustering ... \n')
  clusters_cnts <- as.data.frame(t(as.data.frame(as.vector(table(cluster_vector)))))
  rownames(clusters_cnts) <- NULL
  colnames(clusters_cnts) <- names(table(cluster_vector))
  clusters_cnts <<- clusters_cnts[ , gtools::mixedsort(colnames(clusters_cnts)) ]
  
  # Confusion matrix - cell counts per donor
  cat('Creating confusion matrix for cell counts per donor ... \n')
  cM <- confusionMatrix(paste0(cluster_vector),
                                paste0(sample_vector))
  colnames(cM) <- colnames(cM) %>% str_remove("_ATAC")
  cM <- cM[ gtools::mixedsort(row.names(cM)), ]
  cM <- cM[ , gtools::mixedsort(colnames(cM)) ]
  rownames(cM) <- factor(rownames(cM), levels = rownames(cM))
  colnames(cM) <- factor(colnames(cM), levels = colnames(cM)) # Not working
  
  clust_cM <<- pheatmap::pheatmap(
    mat = cM,
    color = paletteContinuous("whiteBlue"),
    border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
    cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
    treeheight_col = 0,
    treeheight_row = 0,
    angle_col = 45,
    number_color = 'black'
  )
  print(clust_cM)
  
  
  # Create plots for group plot
  cat('\nPlotting UMAPs ... \n\n')
  clusters_UMAP <<- plotEmbedding(ArchRProj = ARCHR_ID, colorBy = "cellColData", 
                                     name = CLUSTERS_ID, 
                                     embedding = UMAP_ID,
                                     plotAs = 'points') + ggtitle('Clusters')
  
  clusters_UMAP_BySample <- plotEmbedding(ArchRProj = ARCHR_ID, colorBy = "cellColData", 
                                          name = "Sample", embedding = UMAP_ID,
                                          plotAs = 'points') + ggtitle(paste0('By Sample'))
  
  
  # Stacked barplots
  cat('\nCreating stacked barplots ... \n')
  cnts_per_donor <- as.data.frame(as.matrix(cM)) %>%
    rownames_to_column("Cluster")
  cnts_per_donor$Cluster <- as.factor(cnts_per_donor$Cluster)
  cnts_per_donor_melt <- reshape2::melt(cnts_per_donor, id = 'Cluster')
  cnts_per_donor_melt$Cluster <- factor(cnts_per_donor_melt$Cluster,
                                        levels = rownames(cM))
  
  # Get the levels for type in the required order - https://stackoverflow.com/questions/22231124
  cnts_per_donor_melt$variable = factor(cnts_per_donor_melt$variable,
                                        levels = SAMPLE_IDs)
  cnts_per_donor_melt = arrange(cnts_per_donor_melt, Cluster, dplyr::desc(variable))
  
  # Calculate percentages
  cnts_per_donor_melt = plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, percent = value/sum(value) * 100)
  
  # Format the labels and calculate their positions
  cnts_per_donor_melt <- plyr::ddply(cnts_per_donor_melt, .(Cluster), transform, pos = (cumsum(value) - 0.5 * value))
  cnts_per_donor_melt$label = paste0(sprintf("%.0f", cnts_per_donor_melt$percent), "%")
  
  # Plot - Note this could also be shown with bars filling plot
  plot_stacked_pct <- ggplot(cnts_per_donor_melt, aes(x = factor(Cluster), y = percent, fill = variable)) +
    geom_bar(position = position_stack(), stat = "identity") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2) +
    theme(legend.position = "none",
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5, angle = 45),
          axis.text.y  = element_text(colour = "#000000", size = 12)) +
    xlab(NULL) + ylab(NULL)
  
  cat('Creating group plot ... \n')
  group_plot <- plot_grid(clusters_UMAP, clusters_UMAP_BySample, plot_stacked_pct,
                          clust_cM$gtable, ncol = 2, align = 'hv', axis = 'rl')
  
  assign(paste0('group_plot_', CLUSTERS_ID), group_plot, .GlobalEnv)
  cat('\nReturned:', paste0('group_plot_', CLUSTERS_ID), '\n\n')
  
  
}

run_initial_qc <- function() {

    # ArchR does some QC when loading the files in so need to load the pre-QC info
    # Pre-filter
    cat('\nLoading pre-filtered data ... \n')
    for (SAMPLE in 1:length(SAMPLE_IDs)) {
    
    # Subset IDs
    sampleID <- substr(SAMPLE_IDs[SAMPLE], 1, 7)
    donorID <- substr(SAMPLE_IDs[SAMPLE], 5, 8)
    
    # Load Pre-filtered data
    preQC_df <- readRDS(paste0(OUT_DIR, "/QualityControl/", SAMPLE_IDs[SAMPLE], "/", 
                                SAMPLE_IDs[SAMPLE], "-Pre-Filter-Metadata.rds"))
    preQC_df$log10nFrags <- log10(preQC_df$nFrags)
    
    # TSS Plot
    cat('\nCreating TSS plot ... \n')
    preQC_tss_uFrag_plot <- ggPoint(
        x = preQC_df[,"log10nFrags"], 
        y = preQC_df[,"TSSEnrichment"], 
        title = SAMPLES[SAMPLE],
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = c(log10(500), quantile(preQC_df[,"log10nFrags"], probs = 0.99)),
        ylim = c(0, quantile(preQC_df[,"TSSEnrichment"], probs = 0.99))
    ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
    
    # Assign frag plots and counts_df
    assign(paste0("preQC_tss_uFrag_plot_", donorID), preQC_tss_uFrag_plot, envir = .GlobalEnv)
    assign(paste0("counts_df_", donorID), 
            data.frame("Sample" = sampleID,
                        "Cells_Pass_Filter" = sum(preQC_df$Keep),
                        "Cells_dropped" = sum(preQC_df$Keep == 0),
                        "Total_Frags" = sum(preQC_df$nFrags),
                        "Median_Frags" = median(preQC_df$nFrags[preQC_df$Keep==1]),
                        "Median_TSS_Enrichment" = median(preQC_df$TSSEnrichment[preQC_df$Keep==1])), envir = .GlobalEnv)
    
    }

    ## Initial QC reporting  --------------------------------------------------------------
    cat('\nCreating pre-filter plots ... \n')

    if (REGION == "FC") {
    
    # Pre-filter tss-frag plot
    cat('\nCreating tss-frag plots ... \n')
    preQC_tss_uFrag_plot <- plot_grid(preQC_tss_uFrag_plot_510, preQC_tss_uFrag_plot_611, preQC_tss_uFrag_plot_993, 
                                        preQC_tss_uFrag_plot_G17, preQC_tss_uFrag_plot_G18, preQC_tss_uFrag_plot_G20, 
                                        preQC_tss_uFrag_plot_G21, preQC_tss_uFrag_plot_T31, preQC_tss_uFrag_plot_T34)
    # Counts df
    cat('\nGenerating counts df ... \n')
    counts_df <<- rbind(counts_df_510, counts_df_611, counts_df_993,
                        counts_df_G17, counts_df_G18, counts_df_G20,
                        counts_df_G21, counts_df_T31, counts_df_T34)
    
    } else {
    
    cat('\nCreating tss-frag plots ... \n')
    preQC_tss_uFrag_plot <- plot_grid(preQC_tss_uFrag_plot_510, preQC_tss_uFrag_plot_611, preQC_tss_uFrag_plot_993)
    # Counts df
    cat('\nGenerating counts df ... \n')
    counts_df <<- rbind(counts_df_510, counts_df_611, counts_df_993)
    
    }

    ## PostQC
    cat('\nCreating post-QC plots ... \n\n')
    archR.meta <- as.data.frame(getCellColData(archR))
    archR.meta$log10nFrags <- log10(archR.meta$nFrags)

    tss_uFrag_plot <<- ggPoint(
    x = archR.meta[,"log10nFrags"], 
    y = archR.meta[,"TSSEnrichment"], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(archR.meta[,"log10nFrags"], probs = 0.99)),
    ylim = c(0, quantile(archR.meta[,"TSSEnrichment"], probs = 0.99))
    ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")


    fragSize_plot <<- plotFragmentSizes(ArchRProj = archR)
    fragSize_plot 

    tss_plot <<- plotTSSEnrichment(ArchRProj = archR)
    tss_plot

    tss_uFrag_plot

    ridge_plot <<- plotGroups(
    ArchRProj =  archR, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
    )

}

# Note this partially contained in run_initial_qc above so maybe redundant code 
print_tss_frag_plot <- function(
  
  ARCHR_ID = NULL ){
    
    #' Creates a tss_frag_plot QC plot for ArchR object 
    #' 
    #' @description Creates a tss_frag_plot QC plot for ArchR object it
    #' requires a variable of SAMPLE_IDs to be set for each sample in the
    #' project
    #' 
    #' @param ARCHR_ID: ArchR object

  
  for (SAMPLE in 1:length(SAMPLE_IDs)) {
    
    cat('\nCreating tss_uFrag_plot for', SAMPLE_IDs[SAMPLE], ' ... \n')
    
    # Subset IDs
    sampleID <- substr(SAMPLE_IDs[SAMPLE], 1, 7)
    donorID <- substr(SAMPLE_IDs[SAMPLE], 5, 8)
    
    QC_df <- as.data.frame(getCellColData(ARCHR_ID)) %>%
      dplyr::filter(Sample == sampleID)
    
    tss_uFrag_plot <- ggPoint(
      x = QC_df[,"log10nFrags"], 
      y = QC_df[,"TSSEnrichment"], 
      title = SAMPLES[SAMPLE],
      colorDensity = TRUE,
      continuousSet = "sambaNight",
      xlabel = "Log10 Unique Fragments",
      ylabel = "TSS Enrichment",
      xlim = c(log10(500), quantile(QC_df[,"log10nFrags"], probs = 0.99)),
      ylim = c(0, quantile(QC_df[,"TSSEnrichment"], probs = 0.99))
    ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
    
    assign(paste0("tss_uFrag_plot_", donorID), tss_uFrag_plot, envir = .GlobalEnv)
    
  }
  
  if (REGION == "FC") {
    
    if (SET_FILTER == 'ziffra_only_default') {
      
      cat('\nReturned tss_uFrag_plot.\n')
      tss_uFrag_plot <<- plot_grid(tss_uFrag_plot_G17, tss_uFrag_plot_G18, tss_uFrag_plot_G20, 
                                   tss_uFrag_plot_G21, tss_uFrag_plot_T31, tss_uFrag_plot_T34)
      
      
    } else {
      
    
    # Pre-filter tss-frag plot
    cat('\nReturned tss_uFrag_plot.\n')
    tss_uFrag_plot <<- plot_grid(tss_uFrag_plot_510, tss_uFrag_plot_611, tss_uFrag_plot_993, 
                                 tss_uFrag_plot_G17, tss_uFrag_plot_G18, tss_uFrag_plot_G20, 
                                 tss_uFrag_plot_G21, tss_uFrag_plot_T31, tss_uFrag_plot_T34)
    
    }
    
  } else {
    
    cat('\nReturned tss_uFrag_plot.\n')
    tss_uFrag_plot <<- plot_grid(tss_uFrag_plot_510, tss_uFrag_plot_611, tss_uFrag_plot_993, 
                                 tss_uFrag_plot_G20, tss_uFrag_plot_T34)
    
    
  }
  
  
  
}

run_cluster_analysis <- function(
    
  ARCHR_ID = NULL,
  LSI_NAME = 'IterativeLSI',
  HARMONY_NAME = NULL,
  COVARIATE = NULL) {
  
  #' Run dimensionality reduction, clustering and generate UMAP 
  #' for archR object and optionally run Harmony batch correction.
  #' 
  #' 
  #' @description Run dimensionality reduction, clustering and generate UMAP 
  #' for archR object and optionally run Harmony batch correction.
  #' Needs to be run seperately for batch correction. Returns archR.2 if
  #' Harmony is not run and archR.3 otherwise.
  #' 
  #' @param ARCHR_ID: ArchR object (required)
  #' @param LSI_NAME: LSI name (required)
  #' @param HARMONY_NAME: Name to use to store Harmony object
  #' @param COVARIATE: The covariate(s) stored in archR object to run batch correction 
  
  if (is.null(HARMONY_NAME)) {
    
    cat('\nRunning dimensionality reduction - pre-batch correction ... \n\n')
    archR.test <- addIterativeLSI(
      ArchRProj = ARCHR_ID,
      useMatrix = "TileMatrix", 
      name = LSI_NAME, 
      iterations = 2, 
      clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10,
        maxClusters = MAX_CLUSTERS
      ), 
      varFeatures = VAR_FEATURES, 
      dimsToUse = 1:30,
      force = TRUE
      
    )
    
    CLUSTERS <- 'Clusters'
    UMAP <- 'UMAP'
    REDUCED_DIM <- LSI_NAME
    
    cat('\nCLUSTERS set to:', CLUSTERS, 
        '\nUMAP set to:', UMAP, 
        '\nREDUCED_DIM set to:', REDUCED_DIM, '... \n\n')
    
  } else {
    
    cat('\nRunning Harmony batch correction ... \n\n')
    archR.test <- addHarmony(
      ArchRProj = ARCHR_ID,
      reducedDims = LSI_NAME,
      name = HARMONY_NAME,
      groupBy = COVARIATE,
      force = TRUE)
    
    CLUSTERS <- 'Clusters_Harmony'
    UMAP <- 'UMAPHarmony'
    REDUCED_DIM <- HARMONY_NAME
    
    cat('\nCLUSTERS set to:', CLUSTERS, 
        '\nUMAP set to:', UMAP, 
        '\nREDUCED_DIM set to:', REDUCED_DIM, 
        '\nCOVARIATE:', COVARIATE, '... \n\n')
    
  }
  
  cat('\nClustering ... \n\n')
  archR.test <- addClusters(
    input = archR.test,
    reducedDims = REDUCED_DIM,
    method = "Seurat",
    name = CLUSTERS,
    resolution = 0.8, 
    force = TRUE
  )
  
  cat('\nCreating UMAP ... \n\n')
  archR.test <- addUMAP(
    ArchRProj = archR.test, 
    reducedDims = REDUCED_DIM, 
    name = UMAP, 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", 
    force = TRUE
  )
  
  if (is.null(HARMONY_NAME)) {
    
    assign('archR.2', archR.test, .GlobalEnv)
    cat('\nReturned archR.2 object.\n')
    
  } else {
    
    assign('archR.3', archR.test, .GlobalEnv)
    cat('\nReturned archR.3 object.\n')
    
    
  }
  
}

run_unconstrained_intergration <- function(
    
  ARCHR_ID = NULL,
  CLUSTERS_ID = NULL,
  UMAP_ID = NULL,
  REDUCED_DIMS = NULL,
  SEURAT_OBJECT = NULL) {
  
  #' Run unconstrained integration on archR object 
  #' 
  #' @description Run unconstrained integration on archR object and Seurat object.
  #' First step of integrating snATACseq and snRNAseq data.
  #' 
  #' @param ARCHR_ID: ArchR object 
  #' @param CLUSTERS_ID: Name of Cluster object stored in ArchR object
  #' @param UMAP_ID: Name of UMAP object stored in ArchR object
  #' @param REDUCED_DIMS: Reduced dims object to run integration on
  #' @param SEURAT_OBJECT: Seurat object 

  
  cat('\nRunning unconstrained integration on', REDUCED_DIMS, '... \n\n')
  archR.test <- addGeneIntegrationMatrix(
    ArchRProj = ARCHR_ID,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = REDUCED_DIMS,
    seRNA = SEURAT_OBJECT,
    addToArrow = FALSE,
    groupRNA = "cellIDs",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un")
  
  cat('\nCreating confusion matrix ... \n')
  cM_geneExp <- as.matrix(confusionMatrix(unname(unlist(getCellColData(archR.test)[CLUSTERS_ID])),
                                          archR.test$predictedGroup_Un)) 
  cM_geneExp <- cM_geneExp[ gtools::mixedsort(row.names(cM_geneExp)), ]
  cM_geneExp <- cM_geneExp[ , gtools::mixedsort(colnames(cM_geneExp)) ]
  rownames(cM_geneExp) <- factor(rownames(cM_geneExp), levels = rownames(cM_geneExp))
  colnames(cM_geneExp) <- factor(colnames(cM_geneExp), levels = colnames(cM_geneExp)) # Not working
  clust_CM_geneExp <- pheatmap::pheatmap(
    mat = as.matrix(cM_geneExp),
    color = paletteContinuous("whiteBlue"),
    border_color = "black", display_numbers = TRUE, number_format =  "%.0f",
    cluster_rows = F, # Needed for row order https://stackoverflow.com/questions/59306714
    treeheight_col = 0,
    treeheight_row = 0,
    angle_col = 45,
    number_color = 'black'
  )
  
  cat('\nCreating cell ID mapping data frame ... \n')
  preClust <- colnames(cM_geneExp)[apply(cM_geneExp, 1 , which.max)]
  integration_df <- t(as.data.frame(cbind(preClust, rownames(cM_geneExp)))) #Assignments
  rownames(integration_df) <- c("RNA", "ATAC")
  colnames(integration_df) <- NULL
  
  cat('\nPlotting UMAPs ... \n')
  atacSeq_UMAP <- plotEmbedding(ArchRProj = archR.test, colorBy = "cellColData", 
                                name = CLUSTERS_ID, embedding = UMAP_ID) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    ggtitle('scATACseq IDs')
  rnaSeq_UMAP <- plotEmbedding(ArchRProj = archR.test, colorBy = "cellColData", 
                               name = "predictedGroup_Un", embedding = UMAP_ID) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    ggtitle('scRNAseq IDs')
  compare_UMAP_plot <- plot_grid(atacSeq_UMAP, rnaSeq_UMAP, align = 'hv')
  
  assign(paste0('clust_CM_geneExp_', CLUSTERS_ID), clust_CM_geneExp, .GlobalEnv)
  assign(paste0('integration_df_', CLUSTERS_ID), integration_df, .GlobalEnv)
  assign(paste0('compare_UMAP_plot_', CLUSTERS_ID), compare_UMAP_plot, .GlobalEnv)
  
  if (REDUCED_DIMS == 'IterativeLSI') {
    
    assign('archR.2', archR.test, .GlobalEnv)
    OBJ_NUM <- '2'
    
  } else {
    
    assign('archR.3', archR.test, .GlobalEnv)
    OBJ_NUM <- '3'
    
  }
  
  cat('\nReturned: archR.', OBJ_NUM, ' object, ',
      paste0('clust_CM_geneExp_', CLUSTERS_ID), ', ',
      paste0('integration_df_', CLUSTERS_ID), ', ',
      paste0('compare_UMAP_plot_', CLUSTERS_ID), '\n', sep = '')
  
}

plot_UMAPs_by_marker_genes <- function( 
    
  ARCHR_ID = NULL,
  UMAP_ID = NULL,
  GENES = NULL) {
  
  #' Plot a set of UMAPs coloured by marker genes
  #' 
  #' @description Plot a set of UMAPs coloured by marker gene
  #' 
  #' @param ARCHR_ID: ArchR object 
  #' @param UMAP_ID: Name of UMAP object stored in ArchR object
  #' @param GENES: A vector of marker genes
  
  archR.test <- addImputeWeights(ARCHR_ID) 
  
  genes_UMAP <- plotEmbedding(
    ArchRProj = archR.test, 
    colorBy = "GeneScoreMatrix", 
    name = GENES, 
    embedding = UMAP_ID,
    imputeWeights = getImputeWeights(archR.test),
    plotAs = 'points',  # See streaking issue https://github.com/GreenleafLab/ArchR/issues/1731
  )
  
  all_genes_UMAP <- lapply(genes_UMAP, function(x){
    
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  
  all_genes_UMAP <- plot_grid(plotlist = all_genes_UMAP)
  
  if (UMAP_ID == 'UMAP') {
    
    assign('all_genes_UMAP', all_genes_UMAP, .GlobalEnv)
    cat('\nReturned: all_genes_UMAP.\n\n')
    
  } else {
    
    assign('all_genes_UMAP_Harmony', all_genes_UMAP, .GlobalEnv)
    cat('\nReturned: all_genes_UMAP_Harmony.\n\n')
    
  }
  
}


# This can't handle Harmony objects yet!!!!
run_peak_calling <- function(
    
  ARCHR_ID = NULL,
  CLUSTERS_ID = NULL,
  NEW_LABELS = NULL,
  OLD_LABELS = NULL,
  MACS2_PATH = NULL,
  FDR_THRESHOLD = 0.05,
  EXTEND_WINDOW = 250) {
  
  # Set macs2 path - note that you need to set the default python env to python 2
  
  cat(paste0('\nAssigning cell IDs to clusters for ', REGION, ' ... \n\n'))
  archR.test <- ARCHR_ID
  CLUSTERS <- unname(unlist(getCellColData(archR.test, CLUSTERS_ID)))
  archR.test$Clusters_broad <- mapLabels(CLUSTERS, 
                                         newLabels = NEW_LABELS, 
                                         oldLabels = OLD_LABELS)
  
  cat(paste0('Create pseudo-bulk replicates for ', REGION, ' ... \n\n'))
  archR.test <- addGroupCoverages(ArchRProj = archR.test, groupBy = "Clusters_broad", force = TRUE)
  
  cat(paste0('\nCalling peaks for ', REGION, ' ... \n'))
  archR.test <- addReproduciblePeakSet(
    ArchRProj = archR.test, 
    groupBy = "Clusters_broad", 
    pathToMacs2 = MACS2_PATH,
    cutOff = FDR_THRESHOLD, 
    extendSummits = EXTEND_WINDOW)
  
  cat(paste0('\nCreate tables and plots for report ', REGION, ' ... \n'))
  coverageParams <- archR.test@projectMetadata$GroupCoverages[["Clusters_broad"]]$Params
  coverage_metadata <- archR.test@projectMetadata$GroupCoverages[["Clusters_broad"]]$coverageMetadata
  maxPeaks_default <- 150000
  peaksPerCell_default <- 500
  
  tableGroups <- table(getCellColData(archR.test, "Clusters_broad", drop = TRUE))
  peakCallParams_summary_df <- lapply(seq_along(coverageParams$cellGroups), function(y){
    x <- coverageParams$cellGroups[[y]]
    uniq <- unique(unlist(x))
    n <- lapply(x, length) %>% unlist %>% sum
    nmin <- lapply(x, length) %>% unlist %>% min
    nmax <- lapply(x, length) %>% unlist %>% max
    data.frame(
      Group=names(coverageParams$cellGroups)[y], 
      nCells=tableGroups[names(coverageParams$cellGroups)[y]], 
      nCellsUsed=length(uniq), 
      nReplicates=length(x), 
      nMin=nmin, 
      nMax=nmax, 
      maxPeaks = min(maxPeaks_default, length(uniq) * peaksPerCell_default)
    )
  }) %>% Reduce("rbind",.)
  
  # Plot peak call summary
  peak_call_summary <- metadata(archR.test@peakSet)$PeakCallSummary
  peak_call_summary_plot <- ggplot(peak_call_summary, 
                                   aes(fill=Var1, y=Freq, x=Group)) + 
    geom_bar(position="stack", stat="identity") +
    viridis::scale_fill_viridis(discrete = T) +
    ggtitle(paste0("Peak call summary for ", REGION, ' at FDR < ',
                   FDR_THRESHOLD, ' with ', EXTEND_WINDOW, 'bp extension')) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill=guide_legend(title="Annotation")) +
    xlab("") +
    ylab(expression("No. of Peaks"~(x10^{"3"})))
  
  assign('archR.4', archR.test, .GlobalEnv)
  assign('peakCallParams_summary_df', peakCallParams_summary_df, .GlobalEnv)
  assign('peak_call_summary_plot', peak_call_summary_plot, .GlobalEnv)
  cat('\nReturned: archR.4, peakCallParams_summary_df, peak_call_summary_plot.\n')
  
}

# Create a venn grob object from Chip peak anno output
pairwise_venn <- function(venn_counts, name1, name2) {
  
  library(ChIPpeakAnno)
  
  #dev.off()
  
  VENN_COUNTS <- venn_counts
  
  left <- unname(VENN_COUNTS$vennCounts[3,3])  + unname(VENN_COUNTS$vennCounts[4,3]) 
  right <- unname(VENN_COUNTS$vennCounts[4,3]) + unname(VENN_COUNTS$vennCounts[2,3])
  union <- unname(VENN_COUNTS$vennCounts[4,3])
  
  VENN <- draw.pairwise.venn(area1 = left, 
                             area2 = right, 
                             cross.area = union,
                             category = c(name1, name2),
                             scaled = FALSE,
                             fill = c("#0073C2FF", "#EFC000FF"),
                             fontfamily = "sans",
                             
                             # Numbers
                             cex = 2,
                             
                             # Set names
                             cat.cex = 1.5,
                             cat.fontface = "bold",
                             cat.default.pos = "outer",
                             cat.pos = c(-10, 10),
                             cat.dist = c(0.055, 0.055),
                             cat.fontfamily = "sans")
  
  return(VENN)
  
  
}

threeway_venn <- function(venn_counts, name1, name2, name3) {
  
  graphics.off()
  
  VENN_COUNTS <- as_tibble(venn_counts$vennCounts[])
  
  left <- sum(VENN_COUNTS[5:7, 5], VENN_COUNTS[8, 4])
  right <- sum(VENN_COUNTS[c(3, 4, 7), 6], VENN_COUNTS[8, 4])
  middle <- sum(VENN_COUNTS[c(2, 4, 6), 7], VENN_COUNTS[8, 4])
  intersect_12 <- sum(VENN_COUNTS[7, 4], VENN_COUNTS[8, 4])
  intersect_23 <- sum(VENN_COUNTS[4, 4], VENN_COUNTS[8, 4])
  intersect_13 <- sum(VENN_COUNTS[6, 4], VENN_COUNTS[8, 4])
  intersect_all3 <- as.double(VENN_COUNTS[8, 4])
  
  VENN <- draw.triple.venn(area1 = left, 
                           area2 = right, 
                           area3 = middle,
                           n12 = intersect_12,
                           n23 = intersect_23, 
                           n13 = intersect_13, 
                           n123 = intersect_all3,
                           category = c(name1, name2, name3),
                           fill = c("#0073C2FF", "#EFC000FF", '#21908DFF'),
                           fontfamily = "sans",
                           
                           # Numbers
                           cex = 2,
                           
                           # Set names
                           cat.cex = 2,
                           cat.fontface = "bold",
                           cat.default.pos = "outer",
                           cat.pos = c(-27, 27, 135),
                           cat.dist = c(0.055, 0.055, 0.085),
                           cat.fontfamily = "sans")
  
  return(VENN)
  
  
}

# To do: function for  LSI param testing

# PARAMETER <- 'Max_Clusters'
# 
# # Iterative LSI parameters
# ITERATIONS <- c(2, 4, 6, 8, 10) # Default 2
# RESOLUTIONS <- c(0.1, 0.2, 0.4, 0.8, 1, 2) # Default 0.2
# VAR_FEATURES <- c(10000, 15000, 20000, 25000, 30000) # Default 25000
# MAX_CLUSTERS <- c(2, 4, 6, 8, 10) # Default 6
# N.STARTS <- c(6, 8, 10, 12) # Default 10
# DIMENSIONS <- c(20, 25, 30, 35) # Default 30
# SAMPLE_CELLS <- c(5000, 7500, 10000) # Default 10000 - crashes when > 10000
# 
# # Set the Variables for LSI
# if (PARAMETER == 'Iteration') {
#   
#   PARAMETERS <- ITERATIONS
#   RESOLUTION <- 0.2
#   VAR_FEATURE <- 25000
#   MAX_CLUSTER <- 6
#   N.START <- 10
#   DIMENSION <- 30
#   SAMPLE_CELL <- 10000
#   
# } else if (PARAMETER == 'Resolution') {
#   
#   ITERATION <- 2
#   PARAMETERS <- RESOLUTIONS
#   VAR_FEATURE <- 25000
#   MAX_CLUSTER <- 6
#   N.START <- 10
#   DIMENSION <- 30
#   SAMPLE_CELL <- 10000
#   
# } else if (PARAMETER == 'Variable_Features') {
#   
#   ITERATION <- 2
#   RESOLUTION <- 0.2
#   PARAMETERS <- VAR_FEATURES 
#   MAX_CLUSTER <- 6
#   N.START <- 10
#   DIMENSION <- 30
#   SAMPLE_CELL <- 10000
#   
# } else if (PARAMETER == 'Max_Clusters') {
#   
#   ITERATION <- 2
#   RESOLUTION <- 0.2
#   VAR_FEATURE <- 25000
#   PARAMETERS <- MAX_CLUSTERS
#   N.START <- 10
#   DIMENSION <- 30
#   SAMPLE_CELL <- 10000
#   
# } else if (PARAMETER == 'N_starts') {
#   
#   ITERATION <- 2
#   RESOLUTION <- 0.2
#   VAR_FEATURE <- 25000
#   MAX_CLUSTER <- 6
#   PARAMETERS <- N.STARTS
#   DIMENSION <- 30
#   SAMPLE_CELL <- 10000
#   
# } else if (PARAMETER == 'Dimensions') {
#   
#   ITERATION <- 2
#   RESOLUTION <- 0.2
#   VAR_FEATURE <- 25000
#   MAX_CLUSTER <- 6
#   N.START <- 10
#   PARAMETERS <- DIMENSIONS
#   SAMPLE_CELL <- 10000
#   
# } else {
#   
#   ITERATION <- 2
#   RESOLUTION <- 0.2
#   VAR_FEATURE <- 25000
#   MAX_CLUSTER <- 6
#   N.START <- 10
#   DIMENSION <- 30
#   PARAMETERS <- SAMPLE_CELLS
#   
# }
# 
# # Loop to set parameter to test 
# for (i in 1:length(PARAMETERS)) {
#   
#   cat(paste0('\nRunning LSI changing ', PARAMETER, ' param to ', PARAMETERS[i], ' ... \n'))
#   
#   if (PARAMETER == 'Iteration') {
#     
#     ITERATION <- PARAMETERS[i]
#     
#   } else if (PARAMETER == 'Resolution') {
#     
#     RESOLUTION <- PARAMETERS[i]
#     
#   } else if (PARAMETER == 'Variable_Features') {
#     
#     VAR_FEATURE <- PARAMETERS[i]
#     
#   } else if (PARAMETER == 'Max_Clusters') {
#     
#     MAX_CLUSTER <- PARAMETERS[i]
#     
#   } else if (PARAMETER == 'N_starts') {
#     
#     N.START <- PARAMETERS[i]
#     
#   } else if (PARAMETER == 'Dimensions') {
#     
#     DIMENSION <- PARAMETERS[i]
#     
#   } else {
#     
#     SAMPLE_CELL <- PARAMETERS[i]
#     
#   }
#   
#   # Set names
#   LSI_NAME <- paste0("LSI_", PARAMETER, "_", PARAMETERS[i])
#   CLUSTERS_NAME <- paste0("Clusters_", PARAMETER, "_", PARAMETERS[i])
#   UMAP_NAME <- paste0("UMAP_", PARAMETER, "_", PARAMETERS[i])
#   
#   cat('Running LSI ... \n')
#   archR.LSI_test <- addIterativeLSI(
#     ArchRProj = archR.2,
#     useMatrix = "TileMatrix", 
#     name = LSI_NAME, 
#     iterations = ITERATION, 
#     clusterParams = list( #See Seurat::FindClusters
#       resolution = RESOLUTION, 
#       sampleCells = SAMPLE_CELL, 
#       n.start = N.START
#     ), 
#     varFeatures = VAR_FEATURE, 
#     dimsToUse = 1:DIMENSION
#   )
#   
#   cat('Assigning clusters ... \n')
#   archR.LSI_test <- addClusters(
#     input = archR.LSI_test,
#     reducedDims = LSI_NAME,
#     method = "Seurat",
#     name = CLUSTERS_NAME,
#     resolution = 0.8,
#     force = TRUE
#   )
#   
#   cat('Creating UMAP ... \n')
#   archR.LSI_test <- addUMAP(
#     ArchRProj = archR.LSI_test, 
#     reducedDims = LSI_NAME, 
#     name = UMAP_NAME, 
#     nNeighbors = 30, 
#     minDist = 0.5, 
#     metric = "cosine"
#   )
#   
#   create_archR_group_plot(archR.LSI_test, CLUSTERS_NAME, UMAP_NAME)
#   create_umap_batch_plots(archR.LSI_test, CLUSTERS_NAME, UMAP_NAME)
#   plot_UMAPs_by_marker_genes(archR.LSI_test, UMAP_NAME, MARKER_GENES)
#   
# }


run_archR2Signac <- function() {
  
  packages <- c("ArchR", "Seurat", "ArchRtoSignac", "Signac", 
                "stringr", "EnsDb.Hsapiens.v86") # required packages
  loadinglibrary(packages)
  
  # Obtain ArchRProject peak matrix for object conversion
  archR <- loadArchRProject(path = paste0(ARCHR_DIR, 'GE/')
  pkm <- getPeakMatrix(archR) # proj is an ArchRProject
  
  # Extract appropriate Ensembl gene annotation and convert to UCSC style.
  # "UCSC" is the default style to change to but can be changed with argument seqStyle
  annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38") 
  
  fragments_dir <- "../results/01CELLRANGER/" # the directory before "/outs/" for all samples
  
  seurat_atac <- ArchR2Signac(
    ArchRProject = archR,
    refversion = "hg38",
    samples = list('14510_WGE_ATAC', '14611_WGE_ATAC', '14993_WGE_ATAC'), 
    fragments_dir = fragments_dir,
    pm = pkm, # peak matrix from getPeakMatrix()
    fragments_fromcellranger = "Yes", # fragments_fromcellranger This is an Yes or No selection ("NO" | "N" | "No" or "YES" | "Y" | "Yes")
    fragments_file_extension = NULL, # Default - NULL: File_Extension for fragments files (typically they should be '.tsv.gz' or '.fragments.tsv.gz')
    annotation = annotations # annotation from getAnnotation()
  )
  
  gsm <- getGeneScoreMatrix(ArchRProject = archR, SeuratObject = seurat_atac)
  
  seurat_atac[['RNA']] <- CreateAssayObject(counts = gsm)
  
  seurat_atac <- addDimRed(
    ArchRProject = archR,
    SeuratObject = seurat_atac,
    addUMAPs = "UMAP",
    reducedDims = "IterativeLSI"
  ) # default is "IterativeLSI"
  
}

