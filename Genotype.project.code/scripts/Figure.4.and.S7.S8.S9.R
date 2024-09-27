### MDM/Monocytes and Neutrophils compartment analysis

# Set project directory
projdir = 'scRNA/MDM.Monocytes/' # define project directory
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
compartment.name = 'MDM.Monocytes'
reductionName = 'sampleID_harmony_umap'


### Proccess MDM.Monocytes compartment
# Load Seurat object
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


saveRDS(srt1, paste0(projDir, 'MDM.Monocytes.srt.rds'))
srt <- srt1
### Reload the compartment if required
srt <- readRDS(paste0(projDir, 'MDM.Monocytes.srt.rds'))
reductionName = 'sampleID_harmony_umap'


##################################################################################
##################################################################################
##################################################################################
### Figure 4a - MDM.Monocytes subcelltype umap. In this figure we are showing the RNA velocity which we calculate the by standerd RNA velocity pipeline by scVelo and we added our embadding with celltype information to show the RNA velocity UMAP with MDM.Monocytes cell susbtype.
metaGroupNames = c('sub_celltype')
umap <- DimPlot (object = srt, reduction = reductionName, pt.size = 0.1, label = FALSE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/sub_celltype_umap_without_label.png'), width = 2200, height = 1500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()


##################################################################################
##################################################################################
##################################################################################
### Figure 4b
meta_col = 'groupID'
output_path = paste0(projDir,'Plots/')
Cell_density_fun(srt= srt, meta_col= meta_col, output_path= output_path, width= 2000, height= 2000)



##################################################################################
##################################################################################
##################################################################################
### Figure 4c - MDM.Monocytes subcelltype cell composition
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


##################################################################################
##################################################################################
##################################################################################
### Figure S7a - Dotplot with canonical marker in MDM/Monocytes
top3.sel.gene <- c( 'Chil3', 'Ccr2', ## Monocytes Marker
                   'Ms4a7', 'Acp5', # MDM marker
    'Cxcl9', 'Cxcl10', 'Rsad2', # Monocytes.Interferon
'Il1b', 'Fn1', 'Tgfbi', # Monocytes.Pro.Inflammatory
'Ahnak', 'S100a4', 'S100a6', # Monocytes.Homeostatic
'Arg1', 'Hilpda', 'Cxcl2', # Monocytes.Disease.Associated
'Gpnmb', 'Spp1', 'Hmox1', # MDM.Disease.Associated
 'Apoe', 'Ctsd', 'Apoc1', # MDM.Homeostatic
'Egr1', 'Fosb', 'Jun', # MDM.Pro.Inflammatory
 'Top2a', 'Stmn1', 'Hist1h1b'# MDM.Proliferating
)
celltype.level <- c('Monocytes/MDM.Interferon', 'Monocytes.Pro.Inflammatory', 'Monocytes.Homeostatic', 
                   'Monocytes.Disease.Associated', 'MDM.Disease.Associated', 'MDM', 'MDM.Pro.Inflammatory', 'MDM.Proliferating')

srt@meta.data <- srt@meta.data |>
     mutate(sub_celltype = factor(sub_celltype, levels = rev(celltype.level)))
p1 <- DotPlot(object = srt, features = top3.sel.gene, scale = T, group.by = 'sub_celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
pdf (paste0(projDir, 'Plots/dotplot_top_3.selected_genes_display_at_sub_celltype.pdf'), useDingbats = F, width = 15, height = 4)
print(p1)
dev.off()


##################################################################################
##################################################################################
##################################################################################
### Figure 4d - wgcna prog boxplot
wgcna_gmt <- GSA.read.gmt(paste0('../../Genotype.prj/data/MDM.Monocytes.wgcna_26modules_25genes_genesets.gmt.txt'))
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

new_colnames = c('Leukocyte.Migration', 'Hypoxia', 'MDM.Cxcl9', 'Extracellular.Space1', 
                 'Positive.Regulation.Of.Cell.Differentiation', 'TNFA.signaling', 
                 'G2M', 'Innate.Immune.Response', 'Unkown1', 'Mono.MDM', 'EMT.KRAS', 'Extracellular.Space2', 
                 'EMT.Myogenesis', 'Interferon', 'TNFA.signaling.Hypoxia', 'Adipogenesis', 'TNFA.signaling', 
                 'Unkown2', 'Ox.Phos', 'E2F.DNA.repair', 'cell.Migration', 
                 'MYC.Target', 'Vacuole', 'KRAS.signaling', 'MYC.Target.MTROC1', 'Homeostatic.Process')

colnames(srt1@meta.data)[(length(colnames(srt@meta.data))+1):length(colnames(srt1@meta.data))] <- new_colnames


srt_wgcna <- srt1
projDirW <- paste0(projDir,'WGCNA.prog/')
dir.create(projDirW)

### Make plot with selected module
new_colnames1 <- c('Hypoxia', 'TNFA.signaling',  
                   'Interferon')
meta_modules_names = new_colnames1
metaGroupNames = c('orig.ident','groupID')
### Generate boxplots per meta groups
ccomp_df = srt_wgcna@meta.data[,new_colnames1, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar
color_bar <- color.list[[metaGroupNames[2]]]
box_p = lapply (seq_along(meta_modules_names), function(x) 
  ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
    geom_boxplot () +
    geom_boxplot(fill = color_bar) +
    geom_jitter (color="black", size=0.4, alpha=0.9) +
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
        plot.title = element_text(face = "bold", color = "black", size = 12)))


pdf (paste0 (projDirW,'WGCNA_module_scores_umap_new_ann_at_groupID_only.boxplot.selected.modules2.pdf'), width = 7, height = 3.3)
print (wrap_plots (box_p, nrow= 1))
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure 4e - Cellphonedb result - Only Spp1 LR selected
## make dotplot with these selected genes
GENE.TO.SHOW <- 'SPP1' # VEGF and TGF
celltype.name <- c('MDM', 'Monocytes')
file.name <- 'MDM.Monocytes'
metaGroupName4 = 'groupID'
genotypes <- unique(srt@meta.data[, metaGroupName4])
genotype_pairs <- combn(genotypes, 2, simplify = FALSE)

df.list <- list()
df.list1 <- list()

for (new_pair_list in genotype_pairs) {
    genotype1 <- new_pair_list[1]
    genotype2 <- new_pair_list[2]
    print(genotype1)
    print(genotype2)
    folder.name <- paste0(genotype1, '_vs_', genotype2, '_cpdb_ds_analysis/')
            input_file_name <- paste0('../../Genotype.prj/data/',folder.name,'lig_rec_differential_highlvl_df.Filtered.Rda')
    df <- get(load(input_file_name))
    pair_compare <- paste0(genotype1, '_vs_', genotype2)


    i <- celltype.name[1]
    df1 <- get(load(input_file_name))
      df1 <- lig_rec_differential_highlvl_df_not_allzeros[
        grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
        grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE), 
      ]

    i <- celltype.name[2]
    df2 <- get(load(input_file_name))
      df2 <- lig_rec_differential_highlvl_df_not_allzeros[
        grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
        grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE), 
      ]
    df <- rbind(df1, df2)




    df$signif <- 0
    df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
    df$signif[df$mut_wt_neglogpval_fishertest > 1.3] <- 2
    df$signif <- factor(df$signif, levels = c(0, 1, 2))

    df1 <- df
    select_best_pvalue <- df[df$mut_wt_neglogpval_fishertest > 0.99,]

    high_freq_neglogpval_fishertest <- names(table(select_best_pvalue$interacting_pair))
    df <- df[df$interacting_pair %in% unique(c(high_freq_neglogpval_fishertest)),]

    df$genotype1 <- genotype1
    df$genotype2 <- genotype2
    df$genotype.pair <- pair_compare
    df.list[[pair_compare]] <- df

    df1$genotype1 <- genotype1
    df1$genotype2 <- genotype2
    df1$genotype.pair <- pair_compare
    df.list1[[pair_compare]] <- df1
}

df_names <- c(names(df.list))
df_names1 <- c(names(df.list1))

comb.df <- Reduce(function(x, y) rbind(x, df.list[[y]]), df_names, init = NULL)
dim(comb.df)
comb.df1 <- Reduce(function(x, y) rbind(x, df.list1[[y]]), df_names1, init = NULL)
dim(comb.df1)
comb.df4 <- comb.df1

comb.df$signif <- as.numeric(as.character(comb.df$signif))
selected_rows <- comb.df[comb.df$signif >= 1,]
comb.df3 <- selected_rows

for (n in 1:dim(comb.df3)[1]) {
    if (comb.df3[n, "mut_wt_prop_diff"] < 0) {
        comb.df3[n, "sig.genotype"] <- comb.df3[n, "genotype2"]
    } else if (comb.df3[n, "mut_wt_prop_diff"] > 0) {
        comb.df3[n, "sig.genotype"] <- comb.df3[n, "genotype1"]
    } else if (comb.df3[n, "mut_wt_prop_diff"] == 0) {
        comb.df3[n, "sig.genotype"] <- 'No.prop.diff'
    }
}
unique.gene.list <- unique(comb.df3$interacting_pair)

filtered_df2 <- data.frame()

for (un.gene.list.name in unique.gene.list) {
    subset_df <- comb.df3[comb.df3$interacting_pair == un.gene.list.name,]
    if (length(unique(subset_df$sig.genotype)) == 1) {
        filtered_df2 <- rbind(filtered_df2, subset_df)
    }
}

interaction.slct <- unique(filtered_df2$interacting_pair)
length(interaction.slct)

un_cell.pair <- unique(filtered_df2$cell.pair)
length(un_cell.pair)

comb.df4 <- comb.df4[comb.df4$cell.pair %in% unique(c(un_cell.pair)),]
comb.df4 <- comb.df4[comb.df4$interacting_pair %in% unique(c(interaction.slct)),]

for (n in 1:dim(comb.df4)[1]) {
    if (comb.df4[n, "mut_wt_prop_diff"] == 0.0) {
        comb.df4[n, "sig.genotype"] <- 'No.prop.diff'
    } else if (comb.df4[n, "mut_wt_prop_diff"] > 0.0) {
        comb.df4[n, "sig.genotype"] <- comb.df4[n, "genotype1"]
    } else if (comb.df4[n, "mut_wt_prop_diff"] < 0.0) {
        comb.df4[n, "sig.genotype"] <- comb.df4[n, "genotype2"]
    }
}

new.list <- unique(comb.df4$sig.genotype)
levels_sequence <- c(new.list[new.list != 'No.prop.diff'], 'No.prop.diff')
comb.df4$sig.genotype <- factor(comb.df4$sig.genotype, levels = levels_sequence)

pdf(paste0(projDir,'Plots/Cellphonedb.dotplot.Spp1.lr.sig.min.2.in.each.row.pdf'), 
height = 1.8, 
width = 4)

comb.df4$signif <- as.numeric(as.character(comb.df4$signif))

egfr.color <- "#0091CA"
nf1.color <- "#D8423D"
pdgfb.color <- "#55AB55"
genotype.color <- c(egfr.color, nf1.color, pdgfb.color)
comb.df4$sig.genotype <- factor(comb.df4$sig.genotype, levels = c('EGFRvIII', 'Nf1', 'PDGFB'))


 gene.name <- GENE.TO.SHOW
# Filter rows where cell.a or cell.b contains any of the selected cell types
comb.df4 <- comb.df4[grepl(gene.name, comb.df4$gene_a) | grepl(gene.name, comb.df4$gene_b), ]

comb.df4 <- comb.df4 %>%
    group_by(interacting_pair, cell.pair) %>%
    filter(sum(signif == 1) + sum(signif == 2) >= 1) %>%
    ungroup()

p <- ggplot(data = comb.df4, aes(x = cell.pair, y = interacting_pair, color = sig.genotype, size = mut_wt_neglogpval_fishertest, fill = sig.genotype)) +
    geom_point(shape = 21, aes(colour = as.factor(signif)), alpha = 0.7) +
    scale_colour_manual(values = c("00FFFFFF", "gray", "black")) +
    scale_fill_manual(values = c(genotype.color, "black")) +
    theme_minimal() +
    facet_wrap(~genotype.pair) +
    labs(
        title = paste0(celltype.name, ' celltype dotplot'),
        x = "Cell Pair",
        y = "Interacting Pair",
        size = "-log10(p-value)",
        color = "Significance",
        fill = "Genotype"
    ) +
    theme(
        text = element_text(size = 5),
        strip.text = element_text(size = 10, face = 'bold'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_size_continuous(range = c(1, 2.5))

print(p)
dev.off()

comb.df4$gene_a <- as.character(comb.df4$gene_a)
comb.df4$gene_b <- as.character(comb.df4$gene_b)

print('Done')



##################################################################################
##################################################################################
##################################################################################
### Figure S8a - Cellphonedb result - Only Spp1 LR selected
## make dotplot with these selected genes
celltype.name <- c('MDM', 'Monocytes')
file.name <- 'MDM.Monocytes'
    metaGroupName4 = 'groupID'
    genotypes <- unique(srt@meta.data[, metaGroupName4])
    genotype_pairs <- combn(genotypes, 2, simplify = FALSE)

    df.list <- list()
    df.list1 <- list()

    for (new_pair_list in genotype_pairs) {
        genotype1 <- new_pair_list[1]
        genotype2 <- new_pair_list[2]
        print(genotype1)
        print(genotype2)
        folder.name <- paste0(genotype1, '_vs_', genotype2, '_cpdb_ds_analysis/')
        input_file_name <- paste0('../../Genotype.prj/data/',folder.name,'lig_rec_differential_highlvl_df.Filtered.Rda')
        df <- get(load(input_file_name))
        pair_compare <- paste0(genotype1, '_vs_', genotype2)

#         i <- celltype.name
#         df <- get(load(input_file_name))
#         df <- lig_rec_differential_highlvl_df_not_allzeros[
#             grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
#             grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE),
#         ]
        i <- celltype.name[1]
        df1 <- get(load(input_file_name))
          df1 <- lig_rec_differential_highlvl_df_not_allzeros[
            grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
            grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE), 
          ]

        i <- celltype.name[2]
        df2 <- get(load(input_file_name))
          df2 <- lig_rec_differential_highlvl_df_not_allzeros[
            grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
            grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE), 
          ]
        df <- rbind(df1, df2)
        
        
        
        
        df$signif <- 0
        df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
        df$signif[df$mut_wt_neglogpval_fishertest > 1.3] <- 2
        df$signif <- factor(df$signif, levels = c(0, 1, 2))

        df1 <- df
        select_best_pvalue <- df[df$mut_wt_neglogpval_fishertest > 0.99,]

        high_freq_neglogpval_fishertest <- names(table(select_best_pvalue$interacting_pair))
        df <- df[df$interacting_pair %in% unique(c(high_freq_neglogpval_fishertest)),]

        df$genotype1 <- genotype1
        df$genotype2 <- genotype2
        df$genotype.pair <- pair_compare
        df.list[[pair_compare]] <- df

        df1$genotype1 <- genotype1
        df1$genotype2 <- genotype2
        df1$genotype.pair <- pair_compare
        df.list1[[pair_compare]] <- df1
    }

    df_names <- c(names(df.list))
    df_names1 <- c(names(df.list1))

    comb.df <- Reduce(function(x, y) rbind(x, df.list[[y]]), df_names, init = NULL)
    dim(comb.df)
    comb.df1 <- Reduce(function(x, y) rbind(x, df.list1[[y]]), df_names1, init = NULL)
    dim(comb.df1)
    comb.df4 <- comb.df1

    comb.df$signif <- as.numeric(as.character(comb.df$signif))
    selected_rows <- comb.df[comb.df$signif >= 1,]
    comb.df3 <- selected_rows

    for (n in 1:dim(comb.df3)[1]) {
        if (comb.df3[n, "mut_wt_prop_diff"] < 0) {
            comb.df3[n, "sig.genotype"] <- comb.df3[n, "genotype2"]
        } else if (comb.df3[n, "mut_wt_prop_diff"] > 0) {
            comb.df3[n, "sig.genotype"] <- comb.df3[n, "genotype1"]
        } else if (comb.df3[n, "mut_wt_prop_diff"] == 0) {
            comb.df3[n, "sig.genotype"] <- 'No.prop.diff'
        }
    }
    unique.gene.list <- unique(comb.df3$interacting_pair)

    filtered_df2 <- data.frame()

    for (un.gene.list.name in unique.gene.list) {
        subset_df <- comb.df3[comb.df3$interacting_pair == un.gene.list.name,]
        if (length(unique(subset_df$sig.genotype)) == 1) {
            filtered_df2 <- rbind(filtered_df2, subset_df)
        }
    }

    interaction.slct <- unique(filtered_df2$interacting_pair)
    length(interaction.slct)

    un_cell.pair <- unique(filtered_df2$cell.pair)
    length(un_cell.pair)

    comb.df4 <- comb.df4[comb.df4$cell.pair %in% unique(c(un_cell.pair)),]
    comb.df4 <- comb.df4[comb.df4$interacting_pair %in% unique(c(interaction.slct)),]

    for (n in 1:dim(comb.df4)[1]) {
        if (comb.df4[n, "mut_wt_prop_diff"] == 0.0) {
            comb.df4[n, "sig.genotype"] <- 'No.prop.diff'
        } else if (comb.df4[n, "mut_wt_prop_diff"] > 0.0) {
            comb.df4[n, "sig.genotype"] <- comb.df4[n, "genotype1"]
        } else if (comb.df4[n, "mut_wt_prop_diff"] < 0.0) {
            comb.df4[n, "sig.genotype"] <- comb.df4[n, "genotype2"]
        }
    }

#     out.fig1 <- paste0(cellphonedb.dir, 'Unique_gene_dotplot_at_celltype.two.significant.level/')
#     dir.create(out.fig1)

    new.list <- unique(comb.df4$sig.genotype)
    levels_sequence <- c(new.list[new.list != 'No.prop.diff'], 'No.prop.diff')
    comb.df4$sig.genotype <- factor(comb.df4$sig.genotype, levels = levels_sequence)

     pdf(paste0(projDir,'Plots/Cellphonedb.dotplot.sig.min.2.in.each.row.pdf'),
    height = length(interaction.slct) / 15,
    width = length(un_cell.pair) / 3)                  
                       

    comb.df4$signif <- as.numeric(as.character(comb.df4$signif))

    egfr.color <- "#0091CA"
    nf1.color <- "#D8423D"
    pdgfb.color <- "#55AB55"
    genotype.color <- c(egfr.color, nf1.color, pdgfb.color)
    comb.df4$sig.genotype <- factor(comb.df4$sig.genotype, levels = c('EGFRvIII', 'Nf1', 'PDGFB'))


    comb.df4 <- comb.df4 %>%
        group_by(interacting_pair, sig.genotype) %>%
        filter(sum(signif >= 1) >= 2) %>%
        ungroup()
                       
    p <- ggplot(data = comb.df4, aes(x = cell.pair, y = interacting_pair, color = sig.genotype, size = mut_wt_neglogpval_fishertest, fill = sig.genotype)) +
        geom_point(shape = 21, aes(colour = as.factor(signif)), alpha = 0.7) +
        scale_colour_manual(values = c("00FFFFFF", "gray", "black")) +
        scale_fill_manual(values = c(genotype.color, "black")) +
        theme_minimal() +
        facet_wrap(~genotype.pair) +
        labs(
            title = paste0(celltype.name, ' celltype dotplot'),
            x = "Cell Pair",
            y = "Interacting Pair",
            size = "-log10(p-value)",
            color = "Significance",
            fill = "Genotype"
        ) +
        theme(
            text = element_text(size = 5),
            strip.text = element_text(size = 10, face = 'bold'),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        ) +
        scale_size_continuous(range = c(1, 2.5))

    print(p)
    dev.off()

    comb.df4$gene_a <- as.character(comb.df4$gene_a)
    comb.df4$gene_b <- as.character(comb.df4$gene_b)
    #write.csv(comb.df4, paste0(out.fig1, file.name, '.df.csv'), row.names = FALSE)
#}
print('Done')






##################################################################################
##################################################################################
##################################################################################
### Figure 4f and S8b - SPP1, CD44, PTGER4 and EGFR expression in TCGA data
gene.name <- c('EGFR')
# Create an empty list
gene_list <- list()
gene_list[['EGFR.gene']] <- gene.name
gene.name <- c('SPP1', 'CD44', 'PTGER4')
for(gene.name in gene.name){
    gene_list[[gene.name]] <- gene.name
}
result <- Add_celltype_expression(gbm.2018, gene_list)
tmp <- result[[1]]
modules <- result[[2]]

result[[1]]$EGFR.gene <- NULL
result[[2]] <- result[[2]][result[[2]] != 'EGFR.gene']


EGFR.color <- "#0091CA"
NF1.color <- "#D8423D"
PDGFB.color <- "#55AB55"
genotype.color <- c(EGFR.color, NF1.color, PDGFB.color)
col.name <- 'new_subtype_based_on_mut_selected_from_supp_7_subtype'
level.factor <- c("EGFRvIII", "NF1", "4q12_PDGFRA")
file.name <- paste0(projDir,"Plots/SPP1.CD44.PTGER4.EGFR.expression_at_NF1.EGFR.4q12_PDGFRA.only.sig.pdf")
width = 4
height = 7
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)


##################################################################################
##################################################################################
##################################################################################
### Figure S8b - EGFR expression in mouse scRNA data -  Figure S8b generated from Tumor cell compartment

























##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
### Neutrophils compartment
### MDM/Monocytes and Neutrophils compartment analysis

# Set project directory
projdir = 'scRNA/Neutrophils/' # define project directory
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
compartment.name = 'Neutrophils'
reductionName = 'sampleID_harmony_umap'

### Proccess Neutrophils compartment
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
saveRDS(srt1, paste0(projDir, 'Neutrophils.srt.rds'))
srt <- srt1
### Reload the compartment if required
srt <- readRDS(paste0(projDir, 'Neutrophils.srt.rds'))
reductionName = 'sampleID_harmony_umap'


##################################################################################
##################################################################################
##################################################################################
### Figure 4g - Neutrophils cell subtype UMAP
metaGroupNames = c('sub_celltype')
umap <- DimPlot (object = srt, reduction = reductionName, pt.size = 0.1, label = FALSE, cols =color.list[[metaGroupNames]], repel = TRUE, group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/sub_celltype_umap.without.label.png'), width = 1600, height = 1200, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()






##################################################################################
##################################################################################
##################################################################################
### Figure 4h - Neutrophils canonical marker dotplot
top.selected.genes <- c('Retnlg', 'Mmp8', 'Lcn2',
                        'Adrb2', 'Hic1', 'Dynll1', 
                        'Fnip2', 'Vegfa', 'Ero1l',
                        'Ifit1', 'Ifit3', 'Cxcl10',
                        'C1qc', 'C1qa', 'Apoe', 
                        'S100a4', 'Ahnak', 'Ccr2',
                        'Rps27l', 'Rpl12', 'Siglecf',
                        'Cxcl2','Cxcl3', 'Ccl6')

cell.lev <- c('Neutrophils.Activated', 'Neutrophils.Adrb2+', 'Neutrophils.Hypoxia', 'Neutrophils.IFN', 'Neutrophils.Macrophage', 'Neutrophils.Monocytes', 'Neutrophils.Ribo', 'Neutrophils.inflammatory')
srt@meta.data <- srt@meta.data |>
     mutate(sub_celltype = factor(sub_celltype, levels = rev(cell.lev)))
p1 <- DotPlot(object = srt, features = top.selected.genes, scale = T, group.by = 'sub_celltype') +
theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
pdf (paste0(projDir, 'Plots/dotplot_top.selected_genes_display_at_sub_celltype.pdf'), useDingbats = F, width = 10, height = 4)
print(p1)
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure 4i - In this figure we are showing the RNA velocity which we calculate the by standerd RNA velocity pipeline by scVelo and we added Microglia compartment embadding with celltype information to show the RNA velocity UMAP.





##################################################################################
##################################################################################
##################################################################################
### Figure S9a - Neutrophils subcelltype cell composition
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


pdf (paste0 (projDir, 'Plots/cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=10, height=5)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()


