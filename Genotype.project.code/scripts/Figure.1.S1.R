
### Main high level compartment analysis

# Set project directory
projdir = 'scRNA/high.level/' # define project directory
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)

# Load Seurat object
srt = readRDS ('../../srt.rds')

## Load TCGA data
gbm.2018 <- readRDS('../../gbm.2018_new.rds')

## Load 
source ('../../PM_scRNA_atlas/scripts/R_libraries.ns.R')
source ('../../PM_scRNA_atlas/scripts/R_utils.ns.R')
source ('../../PM_scRNA_atlas/scripts/palettes.ns.R')

### compartment name
compartment.name = 'High.level'
reductionName = 'sampleID_harmony_umap'

##################################################################################
##################################################################################
##################################################################################
### Figure 1b
celltype.seq <- c('Tumor', 'Stroma', 'Endothelial', 'Astrocyte', 'Oligodendrocytes',  'Excitatory.Neurons', 'Interneurons', 'Choroid.plexus', 'MDM', 'Microglia', 'Monocytes', 'Neutrophils', 'DC', 'pDC',  'Bcells', 'Plasma', 'Tcells', 'NKcells')
srt@meta.data <- srt@meta.data |>
     mutate(celltype = factor(celltype, levels = celltype.seq))

metaGroupNames = c('celltype')
umap <- DimPlot (object = srt, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()


##################################################################################
##################################################################################
##################################################################################
## Figure 1c
selected.gene <- c('Aqp4', 'Slc1a3', 'Cd79a', 'Cd79b', 'Ecrg4', 'Ttr', 'Itgae', 'Xcr1', 'Flt1', 'Vwf', 'Nrgn', 'Camk2a', 'Npy', 'Pmepa1', 'Ms4a7', 'Mrc1', 'P2ry12', 'Tmem119', 'Ccr2', 'Chil3', 'Gzma', 'Klre1', 'S100a9', 'S100a8', 'Apod', 'Mal', 'Jchain', 'Igha', 'Mgp', 'Dcn', 'Cd3d', 'Cd3e', 'Sox11', 'Pdgfra', 'Olig1', 'Olig2', 'Cd300c', 'Ccr9')
celltype.seq <- c('Astrocyte', 'Bcells', 'Choroid.plexus', 'DC', 'Endothelial', 'Excitatory.Neurons', 'Interneurons', 'MDM', 'Microglia', 'Monocytes', 'NKcells', 'Neutrophils', 'Oligodendrocytes', 'Plasma', 'Stroma', 
                 'Tcells', 'Tumor', 'pDC')
srt@meta.data <- srt@meta.data |>
     mutate(celltype = factor(celltype, levels = rev(celltype.seq)))
p1 <- DotPlot(object = srt, features = selected.gene, scale = T, group.by = 'celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) #+
pdf (paste0(projDir,'Plots/Canonical_markers_dotplot.selected.gene1.pdf'), useDingbats = F, width = 13, height = 6)
print(wrap_plots(p1))
dev.off()


##################################################################################
##################################################################################
##################################################################################
## Figure 1d
p1 <- DotPlot(object = srt, features = 'RFP', scale = T, group.by = 'celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral")), position = 'bottom') +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) 

pdf (paste0(projDir,'Plots/RFP.gene_dotplot_at_celltype.pdf'), useDingbats = F, width = 4.5, height = 5)
print(p1)
dev.off()   


##################################################################################
##################################################################################
##################################################################################
## Figure 1e
celltype.seq <- c('Tumor', 'Stroma', 'Endothelial', 'Astrocyte', 'Oligodendrocytes',  'Excitatory.Neurons', 'Interneurons', 'Choroid.plexus', 'MDM', 'Microglia', 'Monocytes', 'Neutrophils', 'DC', 'pDC',  'Bcells', 'Plasma', 'Tcells', 'NKcells')
srt@meta.data <- srt@meta.data |>
     mutate(celltype = factor(celltype, levels = celltype.seq))


metaGroupName1 = 'orig.ident'
metaGroupName2 = 'groupID'
metaGroupName3 = 'celltype'
celltype_length = ceiling(length(unique(srt$celltype))/2)

celltypes_pal_I1 = color.list[[metaGroupName3]]
cc_bar = cellComp (
  seurat_obj = srt,
  metaGroups = c(metaGroupName2,metaGroupName3),
  plot_as = 'bar',
  pal = celltypes_pal_I1
  )

celltypes_pal_I1 = color.list[[metaGroupName2]]
ab <- t(combn(unlist(dimnames(table(srt@meta.data[,metaGroupName2]))),2))
xy.list <- as.list(as.data.frame(t(ab)))
p = cellComp (
  seurat_obj = srt,
  metaGroups = c(metaGroupName1,metaGroupName2,metaGroupName2,metaGroupName3),
  plot_as = 'box',
  pal = celltypes_pal_I1,
  facet_ncol = celltype_length,
  pair_com = xy.list,
  Pvalue_cal = TRUE,
  Pvalue_method = 't_test', # "wilcox_test", "t_test", "sign_test", "dunn_test", "emmeans_test", "tukey_hsd", "games_howell_test"
  hide_ns = TRUE,
  Pvalue_label = "{p.signif}", # "p", "p.adj", "p.signif", "p.adj.signif", "p.format", "p.adj.format"
  p.adjust.method = 'none' # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  )


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=20, height=6)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()



##################################################################################
##################################################################################
##################################################################################
## Figure 1f
## subset only immnune cells
srt1 <- subset(srt, high_level_celltype == 'Lymphoid' | high_level_celltype == 'Myeloid')

## Parameter
celltype = 'only_Immnune_cells'
metaGroupName1 = 'orig.ident'
metaGroupName2 = 'groupID'
metaGroupName3 = 'celltype'
celltype_length = ceiling(length(unique(srt1$celltype))/2)

celltypes_pal_I1 = color.list[[metaGroupName3]]
cc_bar = cellComp (
  seurat_obj = srt1,
  metaGroups = c(metaGroupName2,metaGroupName3),
  plot_as = 'bar',
  pal = celltypes_pal_I1
  )

celltypes_pal_I1 = color.list[[metaGroupName2]]
ab <- t(combn(unlist(dimnames(table(srt1@meta.data[,metaGroupName2]))),2))
xy.list <- as.list(as.data.frame(t(ab)))
p = cellComp (
  seurat_obj = srt1,
  metaGroups = c(metaGroupName1,metaGroupName2,metaGroupName2,metaGroupName3),
  plot_as = 'box',
  pal = celltypes_pal_I1,
  facet_ncol = celltype_length,
  pair_com = xy.list,
  Pvalue_cal = TRUE,
  Pvalue_method = 't_test', # "wilcox_test", "t_test", "sign_test", "dunn_test", "emmeans_test", "tukey_hsd", "games_howell_test"
  hide_ns = TRUE,
  Pvalue_label = "{p.signif}", # "p", "p.adj", "p.signif", "p.adj.signif", "p.format", "p.adj.format"
  p.adjust.method = 'none' # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  )


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj_',celltype,'.pdf'), width=12, height=7)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()



## subset only non-immnune cells
srt1 <- subset(srt, high_level_celltype == 'Astrocyte' | high_level_celltype == 'Endothelial' |
              high_level_celltype == 'Neurons' | high_level_celltype == 'Oligodendrocytes'|
              high_level_celltype == 'Stroma' | high_level_celltype == 'Tumor' | high_level_celltype == 'Choroid.plexus')

## Parameter
celltype = 'only_non_Immnune_cells'
metaGroupName1 = 'orig.ident'
metaGroupName2 = 'groupID'
metaGroupName3 = 'celltype'
celltype_length = ceiling(length(unique(srt1$celltype))/2)

celltypes_pal_I1 = color.list[[metaGroupName3]]
cc_bar = cellComp (
  seurat_obj = srt1,
  metaGroups = c(metaGroupName2,metaGroupName3),
  plot_as = 'bar',
  pal = celltypes_pal_I1
  )

celltypes_pal_I1 = viridis::turbo (nlevels (as.factor (srt1@meta.data[,metaGroupName2])))
celltypes_pal_I1 = color.list[[metaGroupName2]]
ab <- t(combn(unlist(dimnames(table(srt1@meta.data[,metaGroupName2]))),2))
xy.list <- as.list(as.data.frame(t(ab)))
p = cellComp (
  seurat_obj = srt1,
  metaGroups = c(metaGroupName1,metaGroupName2,metaGroupName2,metaGroupName3),
  plot_as = 'box',
  pal = celltypes_pal_I1,
  facet_ncol = celltype_length,
  pair_com = xy.list,
  Pvalue_cal = TRUE,
  Pvalue_method = 't_test', # "wilcox_test", "t_test", "sign_test", "dunn_test", "emmeans_test", "tukey_hsd", "games_howell_test"
  hide_ns = TRUE,
  Pvalue_label = "{p.signif}", # "p", "p.adj", "p.signif", "p.adj.signif", "p.format", "p.adj.format"
  p.adjust.method = 'none' # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  )


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj_',celltype,'.pdf'), width=12, height=7)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()
png (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj_',celltype,'.png'), width=4500, height=2000, res = 300)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure 1g
Neutrophil <- c('CSF3R', 'FPR2', 'IL1R2', 'CXCR2', 'CXCL1')
Microglia <- c('CX3CR1', 'P2RY13', 'P2RY12', 'TMEM119')
Cancer <- c('SOX11', 'ASCL1', 'SOX4', 'HES6', 'DLL3', 'EGFR', 'PTPRZ1', 'FABP7', 'PDGFRA')
# Create a dataframe
gene_data <- data.frame(
  CellType = rep(c('Neutrophil', 'Microglia', 'Cancer'), 
                 times = c(length(Neutrophil), length(Microglia), length(Cancer))),
  Gene = c(Neutrophil, Microglia, Cancer)
)
gene_list <- list()
for(gene_data.name in unique(gene_data$CellType)){
    print(gene_data.name)
    gene_list[[gene_data.name]] <- gene_data$Gene[gene_data$CellType == gene_data.name]
}
result <- Add_celltype_expression(gbm.2018, gene_list)
tmp <- result[[1]]
modules <- result[[2]]
EGFR.color <- "#0091CA"
NF1.color <- "#D8423D"
PDGFB.color <- "#55AB55"
genotype.color <- c(EGFR.color, NF1.color, PDGFB.color)
col.name <- 'new_subtype_based_on_mut_selected_from_supp_7_subtype'
level.factor <- c("EGFRvIII", "NF1", "4q12_PDGFRA")
file.name <- paste0(projDir,"Plots/TCGA.plot_at_NF1.EGFR.4q12_PDGFRA.pdf")
width = 4
height = 8
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)




##################################################################################
##################################################################################
##################################################################################
### Figure S1a
message ('nfeat ncount and mpercentage vln plot')
sampleid_col = 'sampleID'
Idents (srt ) = srt@meta.data[,sampleid_col]
ccomp_df = as.data.frame (table (srt@meta.data[,sampleid_col]))
cc_p2 = ggplot (ccomp_df, aes (x= Var1, y= Freq)) +
    geom_bar (position="stack", stat="identity") +
    theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=0.5, size = 10)) + 
    theme (axis.text = element_text (size = 12)) +
    ggtitle (paste('Total cells',ncol(srt)))+
    theme(plot.title = element_text(hjust = 0.5, size = 15, face="bold"))

vln_p1 = VlnPlot (srt, features = "nFeature_RNA", , group.by = sampleid_col,pt.size = 0, ncol = 1)+
    theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=0.5, size = 10)) +
    theme(legend.position="none")
vln_p2 = VlnPlot (srt, features = "nCount_RNA", group.by = sampleid_col,pt.size = 0, ncol = 1)+
    theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=0.5, size = 10)) +
    theme(legend.position="none")
vln_p3 = VlnPlot (srt, features = "percent.mt", group.by = sampleid_col,pt.size = 0, ncol = 1)+
    theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=0.5, size = 10)) +
    theme(legend.position="none")

png (paste0(projDir,"Plots/QC_nFeat_nCount_m.percent_vlnPlot_at_sample_level.png"), 5000, 1500, res=300)
print ((cc_p2 + vln_p1 + vln_p2 + vln_p3) + plot_layout(widths=c(1,1,1,1)))
dev.off()


##################################################################################
##################################################################################
##################################################################################
### Figure S1b

######################################
## Subset tumor cells
srt1 <- subset(srt, celltype == 'Tumor')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = 'nFeature_RNA'
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                      verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

## Tumor cells UMAP
compartment.name = 'Tumor'
source ('../../Genotype.prj/scripts/palettes.ns.R')
metaGroupNames = 'celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/Tumor.compartment.celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()



######################################
### Microglia compartment
compartment.name = 'Microglia'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'Microglia')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = 'nFeature_RNA'
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

metaGroupNames = 'sub_celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/Microglia.compartment.sub_celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()





######################################
### MDM/Monocytes compartment
compartment.name = 'MDM.Monocytes'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'MDM' | celltype == 'Monocytes')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = 'nFeature_RNA'
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

metaGroupNames = 'sub_celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/MDM.Monocytes.compartment.sub_celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()



######################################
### Neutrophils compartment
compartment.name = 'Neutrophils'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'Neutrophils')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = 'nFeature_RNA'
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

metaGroupNames = 'sub_celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/Neutrophils.compartment.sub_celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()






######################################
### b/plasma compartment
compartment.name = 'Bplasma'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'Bcells' | celltype == 'Plasma')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = NULL
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

metaGroupNames = 'celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/B.Plasma.compartment.celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()


######################################
### TNK compartment
compartment.name = 'TNK'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'Tcells' | celltype == 'NKcells')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = 'nFeature_RNA'
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)


metaGroupNames = 'sub_celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/TNK.compartment.sub_celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()




######################################
### Stroma compartment
compartment.name = 'Stroma'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'Stroma')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = NULL
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

metaGroupNames = 'sub_celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/Stroma.compartment.sub_celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()

######################################
### Endothelial
compartment.name = 'Endothelial'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'Endothelial')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = NULL
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

metaGroupNames = 'sub_celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/Endothelial.compartment.sub_celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()



######################################
### Astro.Ologi.Neuro
compartment.name = 'Astro.Ologi.Neuro'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'Astrocyte' | celltype == 'Choroid.plexus' | celltype == 'Excitatory.Neurons' |                            celltype == 'Interneurons' | celltype == 'Oligodendrocytes')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = NULL
nFeatures = 2000

# Process merged data
srt1 = NormalizeData (object = srt1, normalization.method = "LogNormalize", scale.factor = 10000)
srt1 = FindVariableFeatures (srt1, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1), vars.to.regress = vars_to_regress)
} else {
  srt1 <- ScaleData(srt1, features = VariableFeatures(object = srt1))
}
    
srt1 = RunPCA (srt1, features = VariableFeatures (object = srt1), npcs = ifelse(ncol(srt1) <= 30,ncol(srt1)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt1 = RunUMAP (object = srt1, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt1 = srt1 %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt1 = FindNeighbors (object = srt1, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

metaGroupNames = 'celltype'
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/Astro.Ologi.Neuro.compartment.celltype.png'), width = 2500, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure S1c
### Celldensity UMAP
meta_col = 'groupID'
output_path = paste0(projDir,'Plots/')
Cell_density_fun(srt= srt, meta_col= meta_col, output_path= output_path, width= 2000, height= 2000)