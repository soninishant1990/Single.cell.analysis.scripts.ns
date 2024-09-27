# load libraries
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
  'scWGCNA','igraph','tidyverse', 'paletteer'
)
lapply(packages, require, character.only = TRUE)

library(CellChat)
library(patchwork)
library(cowplot)
library(grid)
library(gridExtra)


# Set cellranger_to_seurat parameters
cr_to_seurat = list(
  org = 'mouse',
  datatype = 'RNA',
  cr_output = 'filtered' # cell range output filter it count matrices
  )

# Set QC parameters (It is cell level filtring)
qc_params = list(
  filtering = 'hard', # 'emptyDrops' or 'hard' filtering
  nFeat = 400, # Number of features per cells. default 400
  nCounts = 1000, # Number of UMI per cell. Default 800
  pchM = 25, # Percent mitochondrial genes. Default 25
  remove.samples = NULL # Remove bad samples. Takes vector of sampleIDs
  )

### Data processing and clustering variables ###
harmony_params = list(
  batch = 'sampleID'
  )

data_processing_param = list(
  nFeatures = 2000, # number of variable genes to consider for dimentionality reduction
  sigPCs = 15, # number of PCA component
  ccRegress = FALSE, # Regress cell cycle gene expression
  metaGroupNames = c('Genotype','sampleID'),
  res = c(0.2, 0.4, 0.8, 2, 3, 5, 8), # denovo cluster resolutions
    vars_to_regress = NULL
  )

# Initiate pipeline
projDir = paste0('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/')
subclusterName='Higher.level.NF1.EGFR.PGDFB'
force = FALSE # re run pipeline from the beginning no matter if objects are found
scrna_pipeline_dir = '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/'
source (paste0(scrna_pipeline_dir,'master_scrna.R'))

srt$high_level_celltype = srt$celltype
srt$high_level_celltype[srt$celltype %in% c('Tumor')] = 'Tumor'
srt$high_level_celltype[srt$celltype %in% c('DC', 'pDC', 'MDM', 'Microglia', 'Monocytes', 'Neutrophils')] = 'Myeloid'
srt$high_level_celltype[srt$celltype %in% c('Bcells', 'NKcells', 'Plasma', 'Tcells')] = 'Lymphoid'
srt$high_level_celltype[srt$celltype %in% c('Stroma')] = 'Stroma'
srt$high_level_celltype[srt$celltype %in% c('Endothelial')] = 'Endothelial'
srt$high_level_celltype[srt$celltype %in% c('Astrocyte')] = 'Astrocyte'
srt$high_level_celltype[srt$celltype %in% c('oligodendrocytes')] = 'Oligodendrocytes'
srt$high_level_celltype[srt$celltype %in% c('Excitatory.Neurons', 'Interneurons')] = 'Neurons'


srt1 <- subset(srt, celltype!='Undefined')
table(srt1$celltype)
srt <- srt1

## Cell chat analysis
## Assign parameter
srt = srt
dir <- paste0(projDir, 'Cellchat_analysis/')
dir.create(dir)
Annotation_level_col <- 'celltype'
high_level_celltype_col = 'high_level_celltype' ## If not put NULL
sample_id_col <- "sampleID"
organism = 'Mouse' # Human

## Assign column for pairwise comparision
Genotype = 'groupID'
Idents(srt) <- Genotype
force = FALSE

## Pathway analysis gene list
pathways.show.list <- c("CXCL")
pathways.show.list <- c("CXCL", "WNT", 'CCL', "TGFb", "VEGF", "SEMA3", "IL1")
#pathways.show.list <- c("CXCL", "WNT", 'CCL', "TGFb", "VEGF", "MHC-I", "MHC-II", 
#"SEMA3", "IFN-I", "SEMA4", "IL1", "IL2")
## get all pathway name - table(CellChatDB$interaction$pathway_name)
# Assuming you have a table named 'interaction_table'
# interaction_table <- table(CellChatDB$interaction$pathway_name)
# # Sort the table in descending order based on the counts
# sorted_table <- sort(interaction_table, decreasing = TRUE)
# # Print the sorted table
# print(sorted_table)


## Create new directory for output file
out_file <- paste0(dir,Annotation_level_col, '/')
dir.create(out_file)
## Create new directory for output file
out_file1 <- paste0(out_file,'Genotype_level_analysis/')
dir.create(out_file1)

## Create new directory for output file
out_file2 <- paste0(out_file,'Genotype_comparison_analysis/')
dir.create(out_file2)


source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/cellchat.fun_single.genotype.analysis.r')
## Run cellchat at genitype level
Cellchat.fun()
## compare at genotype level
genotype_com_fun()


# pathways.show.list <- c("CXCL", "WNT", 'CCL', "TGFb", "VEGF", "SEMA3", "IL1")
# for (Genotype.name in unique(srt@meta.data[, Genotype])){
#   message(Genotype.name)
#   #Genotype.name = 'Nf1'
#   out_file <- paste0(out_file1,Genotype.name, '/')
#   dir.create(out_file)
#   file.name = paste0(out_file, Genotype.name, '.cellchat.object.file.rds')
#   cellchat <- readRDS(file.name)
#   Visualization.of.cell.cell.communication.network()
# }