### Stroma and Endothelial compartment analysis

# Set project directory
projdir = 'scRNA/Stroma/' # define project directory
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
compartment.name = 'Stroma'
reductionName = 'sampleID_harmony_umap'

### Process Stroma compartment
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

saveRDS(srt1, paste0(projDir, 'Stroma.srt.rds'))
srt <- srt1

### Reload the compartment if required
srt <- readRDS(paste0(projDir, 'Stroma.srt.rds'))
reductionName = 'sampleID_harmony_umap'


##################################################################################
##################################################################################
##################################################################################
### Figure 6a - Stroma cell subtype UMAP
metaGroupNames = c('sub_celltype')
umap <- DimPlot (object = srt, reduction = reductionName, pt.size = 1, label = FALSE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/sub_celltype_umap.without.label.png'), width = 2200, height = 1800, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()




##################################################################################
##################################################################################
##################################################################################
### Figure 6b - Stroma subcelltype top 3 marker expression dotplot
cell.seq <- c('Fibroblast', 'Pericytes', 'Smooth.Muscle', 'Stroma.Proliferating')
srt@meta.data <- srt@meta.data |>
    mutate(sub_celltype = factor(sub_celltype, levels = rev(cell.seq)))

selected.marker = c('Dcn','Col1a1','Col1a2',
                   'Higd1b', 'Vtn','Cox4i2',
                    'Acta2','Tagln', 'Mustn1',
                    'Top2a','Ube2c','Mki67')

p1 <- DotPlot(object = srt, features = selected.marker, scale = T, group.by = 'sub_celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) 

pdf (paste0(projDir,'Plots/Stroma_canonical_marker.pdf'), useDingbats = F, width = 8, height = 3.5)
print(p1)
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure 6c-d - Stroma subcelltype cell composition
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


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=8, height=5)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()




##################################################################################
##################################################################################
##################################################################################
### Figure S11  - Daneman et al paper Stroma marker expression in stroma compartment at genotype level
# Pericytes
pericytes <- c('Pdgfrb', 'Cspg4', 'Des', 'Kcnj8', 'Abcc9', 'Anpep')
# Vascular Smooth Muscle
vascular_smooth_muscle <- c('Pdgfrb', 'Acta2', 'Anpep', 'Cspg4', 'Mcam', 'Des')
# Macrophages
macrophages <- c('Cd163', 'Mrc1', 'Lyve1', 'Adgre1')
# Fibroblasts
fibroblasts <- c('Pdgfra', 'Lama1', 'Anpep', 'Col1a1', 'Fn1', 'Pdgfrb')
# Dural Fibroblasts
dural_fibroblasts <- c('Six1', 'Foxp1')
# Arachnoid Fibroblasts
arachnoid_fibroblasts <- c('Crabp2', 'Aldh1a2', 'Slc6a13')
# Pial Fibroblasts
pial_fibroblasts <- c('Ngfr', 'S100a6')

# Combine all markers into a single list
all_markers_caf <- c(pericytes, vascular_smooth_muscle, macrophages, fibroblasts, dural_fibroblasts, arachnoid_fibroblasts, pial_fibroblasts)

# Create a data frame
caf.df <- data.frame(
  celltype = rep(c("Pericytes", "Vascular_Smooth_Muscle", "Macrophages", "Fibroblasts",
                   "Dural_Fibroblasts", "Arachnoid_Fibroblasts", "Pial_Fibroblasts"), 
                 c(length(pericytes), length(vascular_smooth_muscle), length(macrophages), length(fibroblasts), 
                   length(dural_fibroblasts), length(arachnoid_fibroblasts), length(pial_fibroblasts))),
  gene = all_markers_caf
)



srt1 <- srt
for (modules_name in unique(caf.df$celltype)){
    gene_list <- subset(caf.df, celltype == modules_name)
    modulesL = list (gene_list$gene)
    message ('Run AddModuleScore')
    srt1 = AddModuleScore (srt1, modulesL, name = modules_name)
}
meta_modules_names = paste0(unique(caf.df$celltype), '1') # remove grey module (unassigned)
umap_df = data.frame (srt1[[reductionName]]@cell.embeddings, srt1@meta.data[,meta_modules_names])
umap_p1 = lapply (meta_modules_names, function(x) ggplot(data = umap_df) + 
geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = 0.1) + 
scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
 theme_bw() + 
theme(panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank()) +
ggtitle (x))

png (paste0 (projDir,'Plots/Daneman_paper_module_scores_umap.png'), width = 8000, height = 1000, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))))
dev.off()

metaGroupNames = c('orig.ident','groupID')

### Generate boxplots per meta groups
ccomp_df = srt1@meta.data[,meta_modules_names, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt1@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar
color_bar <- color.list[[metaGroupNames[2]]]
box_p = lapply (seq_along(meta_modules_names), function(x) 
  ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
    geom_boxplot () +
    geom_boxplot(fill = color_bar) +
    geom_jitter (color="black", size=0.4, alpha=0.9) +
    geom_pwc(method = "wilcox_test", label = "{p.signif}", hide.ns = T, p.adjust.method = "none") +
    theme_bw() + 
    scale_fill_manual (values= color_bar) + 
    ggtitle (meta_modules_names[x]) + 
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"))+
    theme(plot.title = element_text(size = 8)))

pdf (paste0 (projDir,'Plots/Daneman_paper_module_scores_genes_at_groupID_only.boxplot.pdf'), width = 15, height = 4)
print (wrap_plots (box_p, ncol=7))
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure S11b - Tgfb1 expression in stroma cells
srt1 <- srt
modulesL = list ('Tgfb1')
message ('Run AddModuleScore')
srt1 = AddModuleScore (srt1, modulesL, name = 'Tgfb1')
srt1$Tgfb1 <- srt1$Tgfb11
srt1$Tgfb11 <- NULL
new_colnames1 <- c('Tgfb1')
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

pdf (paste0 (projDir,'Plots/Tgfb1_module_scores.boxplot.pdf'), width = 3, height = 3)
print (wrap_plots (box_p, nrow= 1))
dev.off()




##################################################################################
##################################################################################
##################################################################################
### Figure S11c - Tgfb1 expression in TCGA data
gene.name <- c('TGFB1')
# Create an empty list
gene_list <- list()
for(gene.name in gene.name){
    gene_list[[gene.name]] <- gene.name
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
file.name <- paste0(projDir,"Plots/TGFB1.expression_at_NF1.EGFR.4q12_PDGFRA.only.sig.pdf")
width = 2.5
height = 5
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)








##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
### Endothelial compartment
# Set project directory
projdir = 'scRNA/Endothelial/' # define project directory
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
compartment.name = 'Endothelial'
reductionName = 'sampleID_harmony_umap'


### Process Endothelial compartment
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

saveRDS(srt1, paste0(projDir, 'Endothelial.srt.rds'))
srt <- srt1
### Reload the compartment if required
srt <- readRDS(paste0(projDir, 'Endothelial.srt.rds'))
reductionName = 'sampleID_harmony_umap'


##################################################################################
##################################################################################
##################################################################################
### Figure 6e - Endothelial cell UMAP
metaGroupNames = c('sub_celltype')
umap <- DimPlot (object = srt, reduction = reductionName, pt.size = 1, label = FALSE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/sub_celltype_umap_without_label.png'), width = 2200, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure 6f - Endothelial subcelltype marker dotplot
gene.list <- c('Gm42418', 'AY036118', 'Lars2', 
               'Ackr1', 'Lrg1', 'Cfh', 
               'Cxcl10', 'Ubd', 'Cxcl9',        
               'Cxcl12', 'Slc6a6', 'Mgp', 
               'Plvap', 'Plpp3', 'Plpp1',
               'Slc16a1', 'Spock2', 'Scgb3a1', 
               'Top2a', 'Hmgb2', 'Hist1h2ap', 
               'Angpt2', 'Apln', 'Trp53i11')

celltype.level <- c('EC.Angpt2+', 'EC.CC', 'EC.Cappillary.Venous', 'EC.Choroid.Plexus', 'EC.Cxcl12+', 'EC.Interferon', 'EC.Large.vain', 'EC.Lars2+')

srt@meta.data <- srt@meta.data |>
     mutate(sub_celltype = factor(sub_celltype, levels = celltype.level))

p1 <- DotPlot(object = srt, features = gene.list, scale = T, group.by = 'sub_celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

pdf (paste0(projDir, 'Plots/dotplot_top_3_genes_display_at_sub_celltype.pdf'), useDingbats = F, width = 12, height = 4)
print(p1)
dev.off()




##################################################################################
##################################################################################
##################################################################################
### Figure 6g-h - Endothelial subcelltype cell composition
### Figure 5c-d and S10b - TNK and Bplasma subcelltype cell composition
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


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=12, height=7)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()

