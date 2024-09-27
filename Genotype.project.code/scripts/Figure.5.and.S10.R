### TNK and BPlasma compartment analysis

# Set project directory
projdir = 'scRNA/TNK.bplasma/' # define project directory
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
compartment.name = 'TNK.bplasma'
reductionName = 'sampleID_harmony_umap'

### Process TNK compartment
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

saveRDS(srt1, paste0(projDir, 'TNK.srt.rds'))



##################################################################################
##################################################################################
##################################################################################
### Figure 5a  -TNK cells UMAP
metaGroupNames = c('sub_celltype')
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = FALSE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/TNK.sub_celltype_umap_without_label.png'), width = 2200, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()



##################################################################################
##################################################################################
##################################################################################

compartment.name = 'Bplasma'
source ('../../Genotype.prj/scripts/palettes.ns.R')

srt1 <- subset(srt, celltype == 'Bcells' | celltype == 'Plasma')
srt1$sub_celltype[srt1$sub_celltype == 'B.Plasma.Proliferating'] = 'Plasma'
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

saveRDS(srt1, paste0(projDir, 'Bplasma.srt.rds'))





##################################################################################
##################################################################################
##################################################################################
### Figure S10a  -Bplasma cells UMAP
metaGroupNames = c('sub_celltype')
umap <- DimPlot (object = srt1, reduction = reductionName, pt.size = 0.1, label = FALSE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/Bplasma.sub_celltype_umap_without_label.png'), width = 2200, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()




##################################################################################
##################################################################################
##################################################################################
### Load TNK and Bplasma cells and merge the object for combined analysis
srt.tnk <- readRDS(paste0(projDir, 'TNK.srt.rds'))
srt.tnk
srt.Bplasma <- readRDS(paste0(projDir, 'Bplasma.srt.rds'))
srt.Bplasma

srt.tnk.bplasma <- merge(srt.tnk, y = srt.Bplasma)
srt.tnk.bplasma

saveRDS(srt.tnk.bplasma, paste0(projDir, 'srt.tnk.bplasma.srt.rds'))
srt <- srt.tnk.bplasma


##################################################################################
##################################################################################
##################################################################################
### Figure 5b  -TNK and Bplasma cononical marker dotlpot
desired_order.celltype = c('Bcells', 'Plasma', 'NK', 'NKT','Effector.CD4.T', 'Effector.CD8.T', 'Exhausted.CD4.T', 'Exhausted.CD8.T',
                 'Gamma.delta.T', 'IFN.Responsive.T', 'Memory.T',  'Naive.CD4.T', 'Naive.CD8.T', 'Regulatory.CD4.T', 'T.Ccl4', 'T.Proliferating')

desired_order <- rev(desired_order.celltype)

selected.marker = c('Cd79a','Cd79b',
                    'Jchain', 'Igkc', 
                    'Klre1','Klrk1',
                    'Klra6', 'Cd7',
                    'Cd3d','Cd3e',
                    'Cd4',
                    'Cd8b1','Cd8a',
                    'Gzmk','S100a4',
                    'Tox','Pdcd1',
                    'Trdc', 'Tcrg-C1',
                    'Ifit1','Ifit3',
                     'Lnpep','Mycbp2',
                    'Tcf7','S1pr1',
                    'Foxp3','Ctla4',
                    'Ccl4', 'Nfkb1',
                   'Mki67', 'Top2a')
p1 <- DotPlot(object = srt, features = selected.marker, scale = T, group.by = 'sub_celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
#desired_order <- rev(unique(srt$sub_celltype))
pdf (paste0(projDir, 'Plots/dotplot_canonical_genes_display_at_sub_celltype.pdf'), useDingbats = F, width = 12, height = 4)
print(p1 + scale_y_discrete(limits = desired_order))
dev.off()




##################################################################################
##################################################################################
##################################################################################
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


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=15, height=6)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()



### Figure S10b
srt$sub_celltype1 <- srt$sub_celltype
srt$sub_celltype1[srt$sub_celltype %in% c('Effector.CD8.T', 'Exhausted.CD8.T', 'Naive.CD8.T', 'Effector.CD4.T', 'Exhausted.CD4.T', 'Naive.CD4.T', 'Gamma.delta.T', 'IFN.Responsive.T', 'Memory.T', 'NKT', 'T.Ccl4', 'T.Proliferating', 'Regulatory.CD4.T')] = 'Tcells'
table(srt$sub_celltype1, useNA = "always")

## Parameter
metaGroupName1 = 'orig.ident'
metaGroupName2 = 'groupID'
metaGroupName3 = 'sub_celltype1'
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


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=8, height=3)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()

##################################################################################
##################################################################################
##################################################################################
### Figure 5f - Ctla4 expression in TNK cells
srt1 <- srt.tnk
modulesL = list ('Ctla4')
message ('Run AddModuleScore')
srt1 = AddModuleScore (srt1, modulesL, name = 'Ctla4')

new_colnames1 <- c('Ctla41')
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

pdf (paste0 (projDir,'Plots/Ctla4_module_scores.boxplot.pdf'), width = 3, height = 3)
print (wrap_plots (box_p, nrow= 1))
dev.off()




##################################################################################
##################################################################################
##################################################################################
### Figure 5g - CTLA4 expression in TCGA data
gene.name <- c('CTLA4')
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
file.name <- paste0(projDir,"Plots/CTLA4.expression_at_NF1.EGFR.4q12_PDGFRA.only.sig.pdf")
width = 2.5
height = 4
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)




