#--------------------------------------------------------------------------------------
#
#    ArchR - Initial QC
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Run the analysis up until cluster QC cell removal - Run on Hawk cluster

## Load env variables  ----------------------------------------------------------------
source('scripts/snATACseq_hawk_Renv.R')
source('scripts/snATACseq_functions.R')

##  Load snATACseq data - Cptr 1.5  ---------------------------------------------------
cat('\nCreating Arrow files ... \n')

ArrowFiles <- createArrowFiles(
  inputFiles = c(paste0(DATA_DIR, SAMPLES[1], "/outs/fragments.tsv.gz"),
                 paste0(DATA_DIR, SAMPLES[2], "/outs/fragments.tsv.gz"),
                 paste0(DATA_DIR, SAMPLES[3], "/outs/fragments.tsv.gz")),
  sampleNames = SAMPLE_IDs,
  minTSS = TSS_THRESH, # Dont set this too high because you can always increase later
  minFrags = FRAGS_THRESH, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  QCDir = paste0(OUT_DIR, "/QualityControl"))
  

##  Doublets  - Cptr 2  ---------------------------------------------------------------
cat('\nCalculating Doublet scores ... \n')
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search 
  # with doublet projection.
  LSIMethod = 1,
  outDir = paste0(OUT_DIR, "/QualityControl")
)

##  Create Arrow project  - Cptr 3  ---------------------------------------------------
cat('\nCreate output directory for Arch R project  ... \n')
dir.create(OUT_DIR, recursive = TRUE) # Required archR doesn't create this for you

cat('\nCreating archR project ... \n')
archR <- ArchRProject(ArrowFiles = ArrowFiles, 
                      outputDirectory = OUT_DIR,
                      copyArrows = TRUE # This is recommened so that if you modify 
                      # the Arrow files you have an original copy for later usage.
)

##  Add annotations to objects  -------------------------------------------------------
archR$Donor <- word(archR$Sample, 2, sep = "_")

##  Save and load Arrow project  - Cptr 3.5  ------------------------------------------
cat('\nSaving archR project ... \n')
saveArchRProject(ArchRProj = archR, 
                 outputDirectory = OUT_DIR, 
                 load = FALSE)

# Load project
# archR <- loadArchRProject(path = "")

##  Inital archR QC -------------------------------------------------------------------
run_initial_qc()

## Filter doublets  -------------------------------------------------------------------
cat('\nFiltering doublets ... \n')
archR <- filterDoublets(archR)
doublet_df <- cbind(as.data.frame(table(archR$Sample)), as.data.frame(table(archR$Sample)))
doublet_df[3] <- NULL
doublet_df$cells_removed <- 100 - doublet_df[3] / doublet_df[2] * 100
colnames(doublet_df) <- c("Sample", "Pre_DoubRem", "Post_DoubRem", "pc_cells_removed")
doublet_df

##  Dimensionality reduction  ---------------------------------------------------------
cat('\nRunning dimensionality reduction - pre-batch correction ... \n')
archR <- addIterativeLSI(
  ArchRProj = archR,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10,
    maxClusters = MAX_CLUSTERS
  ), 
  varFeatures = VAR_FEATURES, 
  dimsToUse = 1:30
  
)

##  Clustering  -----------------------------------------------------------------------
cat('\nClustering cells  ... \n')
archR <- addClusters(
  input = archR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

##  Visualisation  --------------------------------------------------------------------
cat('\nCreating UMAP ... \n')
archR <- addUMAP(
  ArchRProj = archR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

## Clustering - reporting  ------------------------------------------------------------
create_archR_group_plot(archR, 'Clusters', 'UMAP')

# Gene specific UMAPs using imputation
# Note that I'm not saving these imputation weights in Save archR section below
# They are for visual cell IDing only at this stage
archR <- addImputeWeights(archR)

genes_UMAP <- plotEmbedding(
  ArchRProj = archR, 
  colorBy = "GeneScoreMatrix", 
  name = MARKER_GENES, 
  embedding = 'UMAP',
  imputeWeights = getImputeWeights(archR),
  plotAs = 'points',  # See streaking issue https://github.com/GreenleafLab/archR/issues/1731
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

save.image(file = paste0('QC_', REGION,'.RData'))

## Save archR project  ----------------------------------------------------------------
cat('\nSaving project ... \n')
saveArchRProject(ArchRProj = archR, 
                 outputDirectory = OUT_DIR, 
                 load = FALSE)

## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
