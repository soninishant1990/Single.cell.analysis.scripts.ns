### Myeloid and Microglia compartment analysis

# Set project directory
projdir = 'scRNA/Myeloid/' # define project directory
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
compartment.name = 'Myeloid'
reductionName = 'sampleID_harmony_umap'


### Proccess myeloid compartment
compartment.name = 'Myeloid'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, high_level_celltype == 'Myeloid')
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


saveRDS(srt1, paste0(projDir, 'Myeloid.srt.rds'))
srt <- srt1
### Reload the compartment if required
srt <- readRDS(paste0(projDir, 'Myeloid.srt.rds'))
reductionName = 'sampleID_harmony_umap'



##################################################################################
##################################################################################
##################################################################################
### Figure 3a - Myeloid cell UMAP
metaGroupNames = c('celltype')
umap <- DimPlot (object = srt, reduction = reductionName, pt.size = 0.1, label = FALSE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/celltype_umap_without_label.png'), width = 2500, height = 1800, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure 3b - cell type cell composition
## Parameter
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


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=9, height=5.5)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()




##################################################################################
##################################################################################
##################################################################################
### Figure S5a - cell type cell composition
srt1 <- srt
cell.seq <- c('pDC', 'DC', 'Neutrophils', 'MDM', 'Monocytes', 'Microglia')
srt1@meta.data <- srt1@meta.data |>
     mutate(celltype = factor(celltype, levels = cell.seq))

selcted.marker.list = c('Sparc', 'P2ry12', 'Tmem119',  
    'Ccr2', 'Chil3', 'F10', 
    'Ms4a7','Mrc1','Acp5',  
    'S100a9','S100a8', 'Cxcr2',
 'Itgae', 'Plet1','Xcr1',
 'Cd300c', 'Cox6a2', 'Klk1')

p1 <- DotPlot(object = srt1, features = selcted.marker.list, scale = T, group.by = 'celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

pdf (paste0(projDir, 'Plots/dotplot_top3.selected.marker.list.pdf'), useDingbats = F, width = 12, height = 4)
print(p1)
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure S5b - cell density plot
meta_col = 'groupID'
output_path = paste0(projDir,'Plots/')
Cell_density_fun(srt= srt, meta_col= meta_col, output_path= output_path, width= 2000, height= 2000)





##################################################################################
##################################################################################
##################################################################################
### Figure S5c - Chmokines 
srt1 = readRDS ('../../srt.rds')
Chemo_gene_list <- c("Ccl2", "Ccl3", "Ccl4", "Ccl5", 'Ccl7', "Ccl8", 'Ccl11',"Ccl12", 
                     'Cxcl1','Cxcl2','Cxcl3','Cxcl5', 'Cxcl9', 'Cxcl10', 'Cxcl11',  
                    "Ccr1", "Ccr2", "Ccr5", "Ccr6", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6",
                    "Cxcr1", "Cxcr2", "Cxcr3", "Cx3cl1")
Chemo_gene_list <- unique(Chemo_gene_list)
srt1 <- subset(srt1, celltype != 'Choroid.plexus')
dotplot_list = list()
for (genoname in names(table(srt1$groupID))){
    srt2 <- subset(srt1, groupID ==genoname)
    p1 <- DotPlot(object = srt2, features = Chemo_gene_list, scale = T, group.by = 'celltype', scale.max = 80) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
    scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
    ggtitle(genoname)
    dotplot_list[[genoname]] <- p1
}

pdf (paste0(projDir, 'Plots/Chemokines_ligand_and_receptor_markers_dotplot_sub_celltype_splited_at_genotype.pdf'), useDingbats = F, width = 30, height = 7)
print(wrap_plots(dotplot_list))
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure 3d - TAMs infiltrating genes expresion
gene.list <- c('Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl7', 'Ccl8', 'Ccl12')
srt1 <- srt
srt1$celltype <- as.character(srt1$celltype)
srt1 <- subset(srt1, celltype == 'Monocytes' | celltype == 'MDM' |celltype == 'Microglia')

srt1 = AddModuleScore (srt1, list(gene.list), name = 'TAM.infiltrate.marker')
meta_modules_names = paste0('TAM.infiltrate.marker', '1') # remove grey module (unassigned)
umap_df = data.frame (srt1[[reductionName]]@cell.embeddings, srt1@meta.data[,meta_modules_names])

metaGroupNames = c('orig.ident','groupID')

### Generate boxplots per meta groups
ccomp_df = srt1@meta.data[,meta_modules_names, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt1@meta.data[,metaGroupNames,drop=F]), mean)

color_bar <- color.list[['groupID']]

box_p = lapply (seq_along(meta_modules_names), function(x) 
  ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
    geom_boxplot(fill = color_bar, outlier.shape = NA) +
    geom_jitter (color="black", size=0.8, alpha=0.9) +
    geom_pwc(method = "wilcox_test", label = "{p.format}{p.signif}", hide.ns = T, p.adjust.method = "none") +
    theme_bw() + 
    scale_fill_manual (values= color_bar) + 
    ggtitle (meta_modules_names[x]) + 
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"))+
    theme(plot.title = element_text(size = 8)))

pdf (paste0 (projDir, 'Plots/Markers.indicating.the.increased.presence.of.TAM_splited_at_genotype_only.boxplot.only.MDM.MONO.MICROGLIA.pdf'), width = 3, height = 3)
print (wrap_plots (box_p, ncol= ifelse (length(box_p) > 4,ceiling(length(box_p)/2),length(box_p))))
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure S5d - TAMs infiltrating genes expresion
gene.list <- c('Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl7', 'Ccl8', 'Ccl12')
srt2 <- srt
for (modules_name in gene.list){
    modulesL = list (modules_name)
    message ('Run AddModuleScore')
    srt2 = AddModuleScore (srt2, modulesL, name = modules_name)
}

# Define old and new column names
old_names <- paste0(gene.list, '1')
new_names <- gene.list
# Get the current column names
current_names <- names(srt2@meta.data)
# Replace old names with new names
current_names[current_names %in% old_names] <- new_names
# Assign the updated column names back to the data frame
names(srt2@meta.data) <- current_names


meta_modules_names = gene.list # remove grey module (unassigned)
umap_df = data.frame (srt2[[reductionName]]@cell.embeddings, srt2@meta.data[,meta_modules_names])

metaGroupNames = c('orig.ident','groupID')

### Generate boxplots per meta groups
ccomp_df = srt2@meta.data[,meta_modules_names, drop=FALSE]
#ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 
ccomp_df = aggregate (ccomp_df, by=as.list(srt2@meta.data[,metaGroupNames,drop=F]), mean)

color_bar <- color.list[['groupID']]

box_p = lapply (seq_along(meta_modules_names), function(x) 
  ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
    #geom_violin (trim=TRUE) +
    #geom_boxplot () +
    geom_boxplot(fill = color_bar, outlier.shape = NA) +
    geom_jitter (color="black", size=0.6, alpha=0.9) +
    geom_pwc(method = "wilcox_test", label = "{p.format}{p.signif}", hide.ns = T, p.adjust.method = "none") +
    theme_bw() + 
    scale_fill_manual (values= color_bar) + 
    ggtitle (meta_modules_names[x]) + 
    #facet_wrap (as.formula(paste("~", metaGroupNames[3]))) +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"))+
    theme(plot.title = element_text(size = 8)))

pdf (paste0 (projDir, 'Plots/Markers.indicating.the.increased.presence.of.TAM_splited_at_genotype_only.boxplot.only.MDM.MONO.MICROGLIA.seperate.pdf'), width = 10, height = 6)
print (wrap_plots (box_p, ncol= ifelse (length(box_p) > 4,ceiling(length(box_p)/2),length(box_p))))
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure 3e - TAMs infiltrating genes expresion in TCGA data
other_gene <- c('CCL2', 'CCL3', 'CCL4', 'CCL5', 'Ccl7', 'CCL8', 'CCL13')

# Create an empty list
gene_list <- list()
gene_list[['Macrophage.infiltrating']] <- other_gene
result <- Add_celltype_expression(gbm.2018, gene_list)
tmp <- result[[1]]
modules <- result[[2]]
EGFR.color <- "#0091CA"
NF1.color <- "#D8423D"
PDGFB.color <- "#55AB55"
genotype.color <- c(EGFR.color, NF1.color, PDGFB.color)
col.name <- 'new_subtype_based_on_mut_selected_from_supp_7_subtype'
level.factor <- c("EGFRvIII", "NF1", "4q12_PDGFRA")
file.name <- paste0(projDir, "Plots/TCGA.DATA.CCL.MACROPHAGE.INFILTRATING.GENE.combined.expression_at_NF1.EGFR.4q12_PDGFRA.pdf")
width = 2
height = 4
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)





##################################################################################
##################################################################################
##################################################################################
### Figure S5d - TAMs infiltrating genes expresion
other_gene <- c('CCL2', 'CCL3', 'CCL4', 'CCL5', 'Ccl7', 'CCL8', 'CCL13')

# Create an empty list
gene_list <- list()
for (gene_name in other_gene){
    gene_list[[gene_name]] <- gene_name
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
file.name <- paste0(projDir, "Plots/TCGA.DATA.CCL.MACROPHAGE.INFILTRATING.GENE.at_NF1.EGFR.4q12_PDGFRA.pdf")
width = 5
height = 9
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)





##################################################################################
##################################################################################
##################################################################################
### Figure 3f - MDM and Monocytes cell's gene expression in TCGA data
library ('clusterProfiler')
gmt.file = paste0 ('/ahg/regevdata/projects/ICA_Lung/Mimi/gbm_atlas/Johnson21/markers.celltypeIMv2GBMF.gmt')
pathways = read.gmt(gmt.file)
pathways <- subset(pathways, term == 'MDM' | term == 'Monocytes')
MDM <- c('APOC1', 'CD163', 'F13A1', 'ITGA4')
Monocytes <- c('S100A9', 'S100A8', 'VCAN', 'FCN1', 'LYZ')
# Create a dataframe
gene_data <- data.frame(
  CellType = rep(c('MDM', 'Monocytes'), 
                 times = c(length(MDM), length(Monocytes))),
  Gene = c(MDM, Monocytes)
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
file.name <- paste0(projDir,"Plots/TCGA.plot_MDM.Monocytes.gene.expression.NF1.EGFR.4q12_PDGFRA.pdf")
width = 4
height = 4
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)


















##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
### Proccess microglia cells
projDir = '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Github.test/scRNA/Microglia/' # define project directory
system (paste('mkdir -p',paste0(projDir,'Plots/')))
setwd (projDir)
# Load Seurat object
srt = readRDS ('../../srt.rds')
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

saveRDS(srt1, paste0(projDir, 'Microglia.srt.rds'))
srt <- srt1
### Load the object if required
srt <- readRDS(paste0(projDir, 'Microglia.srt.rds'))
reductionName = 'sampleID_harmony_umap'


##################################################################################
##################################################################################
##################################################################################
### Figure 3h - Microglia subcelltype canonical markers dotplot
selected.gene.list <- c('Top2a', 'Mki67', 'Stmn1', 
                        'Ifit3', 'Ifit2', 'Ifit1',  
                        'P2ry12', 'Tmem119', 'Cst3',
                        'Gpnmb', 'Spp1', 'Anxa2', 
                       'Il1a', 'Il1b', 'Klf6')
cell.seq <- c('Microglia.Pro.Inflammatory', 'Microglia.Homeostatic', 'Microglia.Disease.Associated', 'Microglia.Interferon.Responsive', 'MIcroglia.Proliferating')

srt@meta.data <- srt@meta.data |>
     mutate(sub_celltype = factor(sub_celltype, levels = cell.seq))


p1 <- DotPlot(object = srt, features = selected.gene.list, scale = T, group.by = 'sub_celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

pdf (paste0(projDir, 'Plots/dotplot_selected.canonical.markers_display_at_sub_celltype.pdf'), useDingbats = F, width = 12, height = 4)
print(p1)
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure 3i - Microglia subcelltype cell composition
## Parameter
metaGroupName1 = 'orig.ident'
metaGroupName2 = 'groupID'
metaGroupName3 = 'sub_celltype'
celltype_length = ceiling(length(unique(srt$sub_celltype))/2)

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


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=13, height=6)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure 3j - In this figure we are showing the RNA velocity which we calculate the by standerd RNA velocity pipeline by scVelo and we added Microglia compartment embadding with celltype information to show the RNA velocity UMAP.


##################################################################################
##################################################################################
##################################################################################
### Figure S6b and S6c - cell density plot
meta_col = 'groupID'
output_path = paste0(projDir,'Plots/')
Cell_density_fun(srt= srt, meta_col= meta_col, output_path= output_path, width= 2000, height= 2000)





##################################################################################
##################################################################################
##################################################################################
### Figure S6b - wgcna prog boxplot
wgcna_gmt <- GSA.read.gmt(paste0('../../Genotype.prj/data/Microglia.wgcna_13modules_25genes_genesets.gmt.txt'))
geneSetNames <- wgcna_gmt$geneset.names
geneSetNames <- gsub("_GeneSet$", "", geneSetNames)
names(wgcna_gmt$genesets) <- geneSetNames
srt1 <- srt
for (modules_name in names(wgcna_gmt$genesets)){
    gene_list <- wgcna_gmt$genesets[[modules_name]][1:25]
    modulesL = list (gene_list)
    message ('Run AddModuleScore')
    srt1 = AddModuleScore (srt1, modulesL, name = modules_name)
}
new_colnames = c('Hypoxia', 
                'Interferon',
                'G2.M',
                'Unknown1', 
                'KRAS.Signaling',
                'XENOBIOTIC.METABOLISM',
                'Allograft.Rejection',
                'TNFA.Signaling', 
                'MTROC1', 
                'Glycolysis',
                'EMT', 
                'Unkown2',
                'Unkown3')

colnames(srt1@meta.data)[(length(colnames(srt@meta.data))+1):length(colnames(srt1@meta.data))] <- new_colnames
metaGroupNames = c('orig.ident','groupID')
projDirW <- paste0(projDir, 'Plots/WGCNA.analysis/')
dir.create(projDirW)

new_colnames1 <- c('Hypoxia', 'Interferon', 'G2.M', 'TNFA.Signaling')
meta_modules_names = new_colnames1
metaGroupNames = c('orig.ident','groupID')
### Generate boxplots per meta groups
ccomp_df = srt1@meta.data[,new_colnames1, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt1@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar
color_bar <- color.list[[metaGroupNames[2]]]
box_p = lapply (seq_along(meta_modules_names), function(x) 
  ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
    geom_boxplot(fill = color_bar, outlier.shape = NA) +
    geom_jitter (color="black", size=0.6, alpha=0.9) +
    geom_pwc(method = "wilcox_test", label = "{p.format}{p.signif}", hide.ns = T, p.adjust.method = "none") +
    theme_bw() + 
    scale_fill_manual (values= color_bar) + 
    ggtitle (meta_modules_names[x]) + 
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", color = "black"),  # Set x-axis text to bold
        axis.text.y = element_text(face = "bold", color = "black"),  # Set y-axis text to bold
        axis.title = element_text(face = "bold", color = "black"),  # Set axis titles to bold
        strip.text = element_text(face = "bold", color = "black"),  # Set facet strip text to bold
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(face = "bold", color = "black", size = 12)
         ))

pdf (paste0 (projDirW,'WGCNA_module_scores_umap_new_ann_at_groupID_only.boxplot.selected.modules.top.25.genes.pdf'), width = 9, height = 3)
print (wrap_plots (box_p, nrow= 1))
dev.off()

srt1$sub_celltype <- as.character(srt1$sub_celltype)
srt1$sub_celltype[srt1$sub_celltype == 'Microglia.Disease.Associated'] = 'MG.DA'
srt1$sub_celltype[srt1$sub_celltype == 'Microglia.Homeostatic'] = 'MG.HOMEO'
srt1$sub_celltype[srt1$sub_celltype == 'Microglia.Interferon.Responsive'] = 'MG.IR'
srt1$sub_celltype[srt1$sub_celltype == 'MIcroglia.Proliferating'] = 'MG.PROLI'
srt1$sub_celltype[srt1$sub_celltype == 'Microglia.Pro.Inflammatory'] = 'MG.PRO.INFLA'
cellseq <- c('MG.PROLI', 'MG.DA', 'MG.HOMEO', 'MG.IR' , 'MG.PRO.INFLA')
srt1@meta.data <- srt1@meta.data |>
     mutate(sub_celltype = factor(sub_celltype, levels = cellseq))

### Make plot with selected module
new_colnames1 <- c('Hypoxia', 'Interferon', 'G2.M', 'TNFA.Signaling')
meta_modules_names = new_colnames1
metaGroupNames = c('orig.ident','sub_celltype')
### Generate boxplots per meta groups
ccomp_df = srt1@meta.data[,new_colnames1, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt1@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar]
color_bar <- color.list[[metaGroupNames[2]]][1:dim(table(srt1@meta.data[,metaGroupNames[2]]))]
box_p = lapply (seq_along(meta_modules_names), function(x) 
  ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
    geom_boxplot(fill = color_bar, outlier.shape = NA) +
    geom_jitter (color="black", size=0.6, alpha=0.9) +
    geom_pwc(method = "wilcox_test", label = "{p.signif}", hide.ns = T, p.adjust.method = "none") +
    theme_bw() + 
    scale_fill_manual (values= color_bar) + 
    ggtitle (meta_modules_names[x]) + 
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", color = "black"),  # Set x-axis text to bold
        axis.text.y = element_text(face = "bold", color = "black"),  # Set y-axis text to bold
        axis.title = element_text(face = "bold", color = "black"),  # Set axis titles to bold
        strip.text = element_text(face = "bold", color = "black"),  # Set facet strip text to bold
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(face = "bold", color = "black", size = 12)
         ))

pdf (paste0 (projDirW,'WGCNA_module_scores_umap_new_ann_at_sub_celltype_only.boxplot.selected.modules.top.25.genes.pdf'), width = 9, height = 3.2)
print (wrap_plots (box_p, nrow= 1))
dev.off()
