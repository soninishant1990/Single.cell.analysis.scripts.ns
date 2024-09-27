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


# srt1 <- subset(srt, celltype!='Undefined')
# table(srt1$celltype)
# srt <- srt1
# table(srt$celltype)
# unique(srt$celltype)

## Cell chat analysis
## Assign parameter
srt = srt
dir <- paste0(projDir, 'Cellchat_analysis/')
dir.create(dir)
Annotation_level_col <- 'celltype'
high_level_celltype_col = 'high_level_celltype' ## If not put NULL
sample_id_col <- "sampleID"
organism = 'Mouse' # Human

## Assign column for pairwise comparision, instead of genotype can provide any other column too
Genotype = 'groupID'
Idents(srt) <- Genotype
force = FALSE

## Pathway analysis gene list
pathways.show.list <- c("CXCL")

## Create new directory for output file
out_file1 <- paste0(dir,Annotation_level_col, '/')
dir.create(out_file1)
## Create new directory for output file
out_file1 <- paste0(out_file1,'Genotype_comparison_analysis/')
dir.create(out_file1)


## merge genotype level cellchat object
out_file2 <- paste0(dir,Annotation_level_col, '/')
dir.create(out_file2)
out_file2 <- paste0(out_file2,'Genotype_level_analysis/')
dir.create(out_file2)

cellchatobj.list <- list()
for (Genotype.name in unique(srt$groupID)){
  message(Genotype.name)
  #Genotype.name = 'Nf1'
  out_file <- paste0(out_file2,Genotype.name, '/')
  file.name = readRDS(paste0(out_file, Genotype.name, '.cellchat.object.file.rds'))
  cellchatobj.list[[Genotype.name]] <- file.name
}

Genotype.name = 'Nf1'
out_file <- paste0(out_file2,Genotype.name, '/')
file.name = readRDS(paste0(out_file, Genotype.name, '.cellchat.object.file.rds'))
file.name1 = subset(file.name, celltype != 'Undefined')



print(cellchatobj.list)
cellchat <- mergeCellChat(cellchatobj.list, add.names = names(cellchatobj.list))


## Save cellchat both object list object and merge object
# Users can now export the merged CellChat object and the list of the two separate objects for later use
save(cellchatobj.list, file = paste0(out_file1, "cellchat_object.list.rds"))
save(cellchat, file = paste0(out_file1, "cellchat_merged.rds"))

out_file <- out_file1
## Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F)
gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight")
png (paste0(out_file,'Compare.the.total.number.of.interactions.and.interaction.strength.png'), width = 1500, height = 1200, pointsize=10, res = 300, type="cairo")
print (wrap_plots (gg1 + gg2))
dev.off()

## Differential number of interactions or interaction strength among different cell populations
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
png (paste0(out_file,'Differential.number.of.interactions.or.interaction.strength.among.different.cell.populations.png'), width = 1500, height = 1200, pointsize=10, res = 300, type="cairo")
print (wrap_plots (gg1 + gg2))
dev.off()

#> Do heatmap based on a merged object
gg1 <- netVisual_heatmap(cellchat, comparison = c(1,2))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
png (paste0(out_file,'Do.heatmap.based.on.a.merged.object.png'), width = 1500, height = 1200, pointsize=10, res = 300, type="cairo")
print (wrap_plots (gg1 + gg2))
dev.off()


source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Higher.level.NF1.EGFR.PGDFB_subset/sampleID_harmony/Cellchat_analysis/celltype/visualization.1.R')

mat1 <-cellchatobj.list[[1]]@net$weight
#mat <- cellchat@net$weight
length(rownames(mat1))

mat2 <-cellchatobj.list[[2]]@net$weight
#mat <- cellchat@net$weight
length(rownames(mat2))

mat3 <-cellchatobj.list[[3]]@net$weight
#mat <- cellchat@net$weight
length(rownames(mat3))


setdiff(union(rownames(mat1), rownames(mat3)), rownames(mat2))











netVisual_circle <- function(net, color.use = NULL, title.name = NULL, sources.use = NULL,
    targets.use = NULL, idents.use = NULL, remove.isolate = FALSE,
    top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL,
    vertex.size.max = NULL, vertex.label.cex = 1, vertex.label.color = "black",
    edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6,
    label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8,
    edge.curved = 0.2, shape = "circle", layout = in_circle(),
    margin = 0.2, vertex.size = NULL, arrow.width = 1, arrow.size = 0.2,
    text.x = 0, text.y = 1.5)
{
    if (!is.null(vertex.size)) {
        warning("'vertex.size' is deprecated. Use `vertex.weight`")
    }
    if (is.null(vertex.size.max)) {
        if (length(unique(vertex.weight)) == 1) {   
            vertex.size.max <- 5
        }
        else {
            vertex.size.max <- 15
        }
    }
    options(warn = -1)
    thresh <- stats::quantile(net, probs = 1 - top) 
    net[net < thresh] <- 0
    if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use))) {
        if (is.null(rownames(net))) {
            stop("The input weighted matrix should have rownames!")
        }
        cells.level <- rownames(net)
        df.net <- reshape2::melt(net, value.name = "value")
        colnames(df.net)[1:2] <- c("source", "target")
        if (!is.null(sources.use)) {
            if (is.numeric(sources.use)) {
                sources.use <- cells.level[sources.use]
            }
            df.net <- subset(df.net, source %in% sources.use)
        }
        if (!is.null(targets.use)) {
            if (is.numeric(targets.use)) {
                targets.use <- cells.level[targets.use]
            }
            df.net <- subset(df.net, target %in% targets.use)
        }
        if (!is.null(idents.use)) {
            if (is.numeric(idents.use)) {
                idents.use <- cells.level[idents.use]
            }
            df.net <- filter(df.net, (source %in% idents.use) |
                (target %in% idents.use))
        }
        df.net$source <- factor(df.net$source, levels = cells.level)
        df.net$target <- factor(df.net$target, levels = cells.level)
        df.net$value[is.na(df.net$value)] <- 0      
        net <- tapply(df.net[["value"]], list(df.net[["source"]],
            df.net[["target"]]), sum)
    }
    net[is.na(net)] <- 0
    if (remove.isolate) {
        idx1 <- which(Matrix::rowSums(net) == 0)    
        idx2 <- which(Matrix::colSums(net) == 0)    
        idx <- intersect(idx1, idx2)
        net <- net[-idx, ]
        net <- net[, -idx]
    }
    g <- graph_from_adjacency_matrix(net, mode = "directed",
        weighted = T)
    edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
    coords <- layout_(g, layout)
    if (nrow(coords) != 1) {
        coords_scale = scale(coords)
    }
    else {
        coords_scale <- coords
    }
    if (is.null(color.use)) {
        color.use = scPalette(length(igraph::V(g))) 
    }
    if (is.null(vertex.weight.max)) {
        vertex.weight.max <- max(vertex.weight)     
    }
    vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
        5
    loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g),
        2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g),
        2]/coords_scale[igraph::V(g), 1]))
    igraph::V(g)$size <- vertex.weight
    igraph::V(g)$color <- color.use[igraph::V(g)]   
    igraph::V(g)$frame.color <- color.use[igraph::V(g)]
    igraph::V(g)$label.color <- vertex.label.color  
    igraph::V(g)$label.cex <- vertex.label.cex      
    if (label.edge) {
        igraph::E(g)$label <- igraph::E(g)$weight   
        igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
    }
    if (is.null(edge.weight.max)) {
        edge.weight.max <- max(igraph::E(g)$weight) 
    }
    if (weight.scale == TRUE) {
        igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max *
            edge.width.max
    }
    else {
        igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
    }
    igraph::E(g)$arrow.width <- arrow.width
    igraph::E(g)$arrow.size <- arrow.size
    igraph::E(g)$label.color <- edge.label.color
    igraph::E(g)$label.cex <- edge.label.cex
    igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,
        1]], alpha.edge)
    igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
    if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
        igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[,
            1])] <- loop.angle[edge.start[which(edge.start[,
            2] == edge.start[, 1]), 1]]
    }
    radian.rescale <- function(x, start = 0, direction = 1) {
        c.rotate <- function(x) (x + start)%%(2 * pi) * direction
        c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
    }
    label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
        direction = -1, start = 0)
    label.dist <- vertex.weight/max(vertex.weight) + 2
    plot(g, edge.curved = edge.curved, vertex.shape = shape,
        layout = coords_scale, margin = margin, vertex.label.dist = label.dist,
        vertex.label.degree = label.locs, vertex.label.family = "Helvetica",
        edge.label.family = "Helvetica")
    if (!is.null(title.name)) {
        text(text.x, text.y, title.name, cex = 1.1)
    }
    gg <- recordPlot()
    return(gg)
}
