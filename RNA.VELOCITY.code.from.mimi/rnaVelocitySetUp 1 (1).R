message ('Load R packages')
packages = c('ggplot2','Seurat','dplyr','Matrix','ggpubr',
  'biomaRt','gplots','patchwork','ComplexHeatmap',
  'RColorBrewer','ggrepel','fgsea',
  'tidyr','gdata','GSA',
  'DropletUtils', 
  'harmony','scater',
  'destiny','plotly',
  'GO.db', 'org.Mm.eg.db','org.Hs.eg.db',
  'viridis',
  'BiocParallel',
  'ggridges',
  'clusterProfiler',
  'NMF',
  'WGCNA',
  'hdWGCNA','igraph','tidyverse', 'hdf5r'
)
lapply(packages, require, character.only = TRUE)

## Set Up Variables
# Where this code is and the metadata will be stores:
proj.dir = 

# Seurat object with embedding saved
srt <- readRDS("srt.rds")
# Name for the metadata file you are creating here
outputFile <- ""

# Add embedding data to metadata: can change umap to harmony variable if needed
srt$umap_1 <- srt@reductions$umap@cell.embeddings[,1]
srt$umap_2 <- srt@reductions$umap@cell.embeddings[,2]

## Creating a barcode variable to match the barcode rownames in the fastq/loom files
srt$barcode <- colnames(srt)

# Seperating out the meta.data
md = srt@meta.data

## INSERT CODE HERE TO EDIT ROWNAMES FROM SEURAT TO MATCH BARCODES IN LOOM FILE
# For example: these are the changes I made for Romain's data originally so it might also be what you need
# md$barcode = gsub("-1", "", md$barcode) # Removes -1 from rownames so they match barcodes
# md$barcode = gsub(".*\\_", "",  md$barcode) # Removes everything before the underscore so rownames match the barcode
# md$barcode = paste0("AMLU01_", md$mouseID, "_0_", md$barcode) # Adds in more information to match barcodes

# Putting meta.data back
srt@meta.data = md

## Can check that it worked
# unique(srt$sampleID)
# tail(srt@meta.data)

# save embedding metadata table
write.csv(srt@meta.data, file= paste0(proj.dir, outputFile), quote=F, row.names=F)