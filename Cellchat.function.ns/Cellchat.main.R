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

###################
### IMPORTANT NOTE - All genotype should have same length of celltype, otherwise it will throw error
################
## Cell chat analysis
## Assign parameter
projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Higher.level.NF1.EGFR.PGDFB_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir, 'srt.rds'))
srt = srt
projDir = projDir
dir <- paste0(projDir, 'Cellchat_analysis_new/') ## create new directory for cellchat
dir.create(dir)
Annotation_level_col <- 'celltype' ## celltype column name
high_level_celltype_col = 'high_level_celltype' ## If not use NULL
sample_id_col <- "sampleID" ## sample id column name
organism = 'Mouse' # Human or Mouse

## Assign column for pairwise comparison
Genotype = 'groupID' ### Here, we can provide the genotype and condition column name
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

## Create new directory for the output file
Geno_file <- paste0(dir,Annotation_level_col, '/')
dir.create(Geno_file)
## Create new directory for the output file
out_file1 <- paste0(Geno_file,'Genotype_level_analysis/')
dir.create(out_file1)

## Create new directory for the output file
out_file2 <- paste0(Geno_file,'Genotype_comparison_analysis/')
dir.create(out_file2)

source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/cellchat.fun.r')
## Run cellchat at genitype level
Cellchat.fun(srt = srt, Genotype = Genotype, out_file1=out_file1, force = force, 
organism = organism, pathways.show.list = pathways.show.list)
## compare at genotype level
genotype_com_fun(out_file1 = out_file1, out_file2 = out_file2, 
Genotype = Genotype, force = force, pathways.show.list = pathways.show.list)

## Make a single file for chord diagram for all pairwise comparisons. 
## Before generating chord diagram it controls the edge weights across different datasets
Comparison_analysis_of_multiple_datasets(out_file1 = out_file1, out_file2 = out_file2, 
Genotype = Genotype, force = force, pathways.show.list = pathways.show.list)

# ### Make pathway cord diagram at 
# file.name = paste0(out_file1, "cellchat_object.list.rds")
# load(paste0(out_file1, "cellchat_object.list.rds")) ## load by cellchatobj.list
# load(paste0(out_file1, "cellchat_merged.rds")) ## load by cellchat

# nf1.obj <- cellchatobj.list[['Nf1']]
# nf1.obj <- updateCellChat(nf1.obj)
# levels(nf1.obj@idents)

# pdgfb.obj <- cellchatobj.list[['PDGFB']]
# pdgfb.obj <- updateCellChat(pdgfb.obj)
# levels(pdgfb.obj@idents)

# egfr.obj <- cellchatobj.list[['EGFRvIII']]
# egfr.obj <- updateCellChat(egfr.obj)
# levels(egfr.obj@idents)

# cellchat <- mergeCellChat(cellchatobj.list, add.names = names(cellchatobj.list), cell.prefix = TRUE)

# out_file5 <- paste0(Geno_file,'Comparison.analysis.of.multiple.datasets.with.different.cell.type.compositions/')
# dir.create(out_file5)
# # Hierarchy plot
# ## get the all pathway list in cellchatDB
# # interaction_table <- table(CellChatDB$interaction$pathway_name)
# pathways.show.list <- c("CXCL", "WNT", 'CCL', "TGFb", "VEGF", "MHC-I", "MHC-II", 
# "SEMA3", "IFN-I", "SEMA4", "IL1", "IL2")
# for (pathways.show in  pathways.show.list){
#   tryCatch({
#     weight.max <- getMaxWeight(cellchatobj.list, slot.name = c("netP"), attribute = pathways.show)
#     plot.list <- list()
#     vertex.receiver = seq(1,10)
#     ## Visualize the inferred signaling network
#     ### Hierarchy plot
#     pdf(paste0(out_file5, pathways.show, '_Visualize.the.inferred.signaling.network.at_genotype_level.Hierarchy.plot.pdf'), width = 8, height = 8)
#     for (i in 1:length(cellchatobj.list)) {
#       p1 <- netVisual_aggregate(cellchatobj.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchatobj.list)[i]))
#       grid.text(paste0("Number of interactions - ", pathways.show, " - ",names(cellchatobj.list)[[i]]), x = 0.2, y = 0.88, gp = gpar(fontsize = 7, fontface = "bold"))
#       print(p1)
#     }
#     dev.off()

#     # Circle plot
#     pdf(paste0(out_file5, pathways.show, '_Visualize.the.inferred.signaling.network.at_genotype_level.Circle.plot.pdf'), width = 8, height = 8)
#     for (i in 1:length(cellchatobj.list)) {
#       p1 <- netVisual_aggregate(cellchatobj.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchatobj.list)[i]))
#       grid.text(paste0("Number of interactions - ", pathways.show, " - ",names(cellchatobj.list)[[i]]), x = 0.2, y = 0.88, gp = gpar(fontsize = 7, fontface = "bold"))
#       print(p1)
#     }
#     dev.off()

#     # Chord diagram
#     pdf(paste0(out_file5, pathways.show, '_Visualize.the.inferred.signaling.network.at_genotype_level.Chord.diagram.pdf'), width = 8, height = 8)
#     for (i in 1:length(cellchatobj.list)) {
#       p1 <- netVisual_aggregate(cellchatobj.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchatobj.list)[i]))
#       grid.text(paste0("Number of interactions - ", pathways.show, " - ",names(cellchatobj.list)[[i]]), x = 0.2, y = 0.88, gp = gpar(fontsize = 7, fontface = "bold"))
#       print(p1)
#     }
#     dev.off()

#   }, error = function(e) {
#       # Handle the error gracefully
#       cat("Error occurred for source", pathways.show, ":", conditionMessage(e), "\n")
#       # You can choose to skip this iteration or take other appropriate actions
#   })
# }











srt<-readRDS("/ahg/regevdata/projects/ICA_Lung/Komal/Lillian/srt_merged_allcompartments.rds")
## Assign parameter
srt = srt
projDir = "/ahg/regevdata/projects/ICA_Lung/Komal/Lillian/"
dir <- paste0(projDir, 'Cellchat_analysis.NS/')
dir.create(dir)
Annotation_level_col <- 'celltype_CPDB'
high_level_celltype_col = 'celltype' ## If not put NULL
sample_id_col <- "sampleID"
organism = 'Mouse' # Human

## Assign column for pairwise comparision
Genotype = 'condition'
Idents(srt) <- Genotype
force = FALSE

## Pathway analysis gene list
pathways.show.list <- c("CXCL", "WNT", 'CCL', "TGFb", "VEGF", "SEMA3", "IL1")

## Create new directory for output file
Geno_file <- paste0(dir,Annotation_level_col, '/')
dir.create(Geno_file)
## Create new directory for output file
out_file1 <- paste0(Geno_file,'Genotype_level_analysis/')
dir.create(out_file1)

## Create new directory for output file
out_file2 <- paste0(Geno_file,'Genotype_comparison_analysis/')
dir.create(out_file2)


source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/cellchat.fun.r')
## Run cellchat at genitype level
Cellchat.fun(srt = srt, Genotype = Genotype, out_file1=out_file1, force = force, 
organism = organism, pathways.show.list = pathways.show.list)
## compare at genotype level
genotype_com_fun(out_file1 = out_file1, out_file2 = out_file2, 
Genotype = Genotype, force = force, pathways.show.list = pathways.show.list)

## Make a single file for chord diagram for all paiwise comparision. Before generating chord 
## diagram it control the edge weights across different datasets
Comparison_analysis_of_multiple_datasets(out_file1 = out_file1, out_file2 = out_file2, 
Genotype = Genotype, force = force, pathways.show.list = pathways.show.list)
