### Tumor compartment analysis

# Set project directory
projdir = 'scRNA/Tumor/' # define project directory
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
compartment.name = 'Tumor'
reductionName = 'sampleID_harmony_umap'


### Proccess Tumor compartment
srt <- subset(srt, celltype == 'Tumor')
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraph = paste0 (paste(batch,collapse='_'),'_harmony_snn')

sigPCs = 15
vars_to_regress = 'nFeature_RNA'
nFeatures = 2000

# Process merged data
srt = NormalizeData (object = srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt = FindVariableFeatures (srt, selection.method = "vst", nfeatures = nFeatures)
if (!is.null(vars_to_regress)) {
  srt <- ScaleData(srt, features = VariableFeatures(object = srt), vars.to.regress = vars_to_regress)
} else {
  srt <- ScaleData(srt, features = VariableFeatures(object = srt))
}
    
srt = RunPCA (srt, features = VariableFeatures (object = srt), npcs = ifelse(ncol(srt) <= 30,ncol(srt)-1,30), ndims.print = 1:5, nfeatures.print = 5, verbose = FALSE)
  
if (batch == 'no')
    {
    srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)
    } else {
    # Run Harmony
    srt = srt %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }

# Run denovo clustering on non-adjusted reductions
srt = FindNeighbors (object = srt, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=reductionGraph)

saveRDS(srt, paste0(projDir, 'Tumor.srt.rds'))

## Load if tumor seurat object is not loaded
srt <- readRDS(paste0(projDir, 'Tumor.srt.rds'))
reductionName = 'sampleID_harmony_umap'


##################################################################################
##################################################################################
##################################################################################
### Figure 1a - cell cycle classify UMAP
message ('Cell cycle analysis')
if (org == 'mouse') {
    m.g2m.genes = readRDS ('../../Genotype.prj/data/cellcycle_mouse.g2m.genes.rds')
    m.s.genes = readRDS ('../../Genotype.prj/data/cellcycle_mouse.s.genes.rds')
    s.genes <- m.s.genes
    g2m.genes <- m.g2m.genes
} else if (org == 'human'){
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
}
## Calculate CellCycleScoring score
srt <- CellCycleScoring(srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

 ### Thrushold value
ccs_threshold = 0.01
ccg2m_threshold = 0.01

## Classify the cells into cycling and non-cycling
srt$Phase2 = ifelse (srt$S.Score > ccs_threshold | srt$G2M.Score > ccg2m_threshold, 'cycling', 'non-cycle')

### Create UMAP
## Cycling color
non_cycling_color <- "#CCCCCC"  # Light Gray for non-cycling cells
cycling_color <- "#66c2a5"  # Green for cycling cells
cc.color <- c(cycling_color, non_cycling_color)
metaGroupNames = c('Phase2')
umap <- DimPlot (object = srt, reduction = reductionName, pt.size = .1, cols =cc.color, group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/classify_cell_cycle_umap.png'), width = 2000, height = 1800, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()





##################################################################################
##################################################################################
##################################################################################
### Figure S2a - cell cycle classify UMAP
metaGroupNames = c('groupID')
umap <- DimPlot (object = srt, reduction = reductionName, pt.size = .01, label = FALSE, repel = TRUE, cols =color.list[[metaGroupNames]], group.by = metaGroupNames) #+theme(legend.position="bottom")
png (paste0(projDir,'Plots/groupID_umap.without.label.png'), width = 2800, height = 2000, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap))
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure 1a - Cell density umaps
### Celldensity UMAP
meta_col = 'groupID'
output_path = paste0(projDir,'Plots/')
Cell_density_fun(srt= srt, meta_col= meta_col, output_path= output_path, width= 2000, height= 2000)


##################################################################################
##################################################################################
##################################################################################
### Figure 1b - Cell cycle cell composition
## Parameter
metaGroupName1 = 'orig.ident'
metaGroupName2 = 'groupID'
metaGroupName3 = 'Phase2'
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


pdf (paste0 (projDir, 'Plots/cell_composition_cycling_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=6, height=5.5)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()


##################################################################################
##################################################################################
##################################################################################
### Figure 1c - Cell cycle cell composition - TCGA data
Cycling <- c('TOP2A', 'MKI67', 'UBE2C', 'PTTG1', 'HMGB2')
# Create a dataframe
gene_data <- data.frame(
  CellType = rep(c('Cycling'), 
                 times = c(length(Cycling))),
  Gene = c(Cycling)
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
file.name <- paste0(projDir,"Plots/TCGA.plot_Cycling.gene.at_NF1.EGFR.4q12_PDGFRA.pdf")
width = 3
height = 4
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)


##################################################################################
##################################################################################
##################################################################################
### Figure 1e and S2c - Cellphonedb analysis
# Define the genotypes for comparison
genotypes <- unique(srt@meta.data[,Genotype])

# Create a pairwise comparison of genotypes
pair_list <- combn(genotypes, 2, simplify = FALSE)

for (pair in pair_list) {
    genotype1 <- pair[1]
    genotype2 <- pair[2]
    new_pair_list <- c(genotype1, genotype2)

    n = 0
    for (group_name in new_pair_list){
        if (n==0){
            sub_srt <- subset(srt, idents = group_name)
        } else{
            sub_srt1 <- subset(srt, idents = group_name)
            sub_srt <- merge(sub_srt, sub_srt1)
        }
        n = n+1
    }
    print(table(sub_srt@meta.data[,Genotype]))
    sub_srt$genotype_group = ifelse (sub_srt@meta.data[,Genotype] == genotype1,genotype1,genotype2)
    print(table(sub_srt$genotype_group))

    ## Provide the directory where cellphonedb results are store for downstaream analysis
    input.dir <- paste0('Cellphonedb.output.result/')
    sample.names <- list.files(paste0('Cellphonedb.output.result/'))
    sample_col_name <- sample_id_col
    genotype_colname <- Genotype ## 
    genotype_analysis <- c(genotype1,genotype2) ## If
    sub_srt <- sub_srt

    out.dir <- paste0('../../Genotype.prj/data/',genotype1,'_vs_',genotype2,'_cpdb_ds_analysis/')
    dir.create(out.dir)

    ### Generate interaction file for generating dotplot
    Cellphonedb_ds_analysis(input.dir=input.dir, out.dir = out.dir, sample.names = sample.names, 
                                        sample_col_name = sample_col_name, genotype_colname =genotype_colname,
                                        genotype_analysis = genotype_analysis)
    message(new_pair_list)
}



##################################################################################
##################################################################################
##################################################################################
## Figure S2c
### Make dotplot with cellphonedb result
celltype.name = 'Tumor'
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

    i <- celltype.name
    df <- get(load(input_file_name))
    df <- lig_rec_differential_highlvl_df_not_allzeros[
        grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
        grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE),
    ]
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
write.csv(comb.df4, paste0(projDir,'Plots/', celltype.name, '.df.csv'), row.names = FALSE)


### Plot only with WNT ligand receptor
gene.name <- 'WNT'
# Filter rows where cell.a or cell.b contains any of the selected cell types
comb.df4 <- comb.df4[grepl(gene.name, comb.df4$gene_a) | grepl(gene.name, comb.df4$gene_b), ]
comb.df4$sig.genotype = factor(comb.df4$sig.genotype, levels =c('EGFRvIII', 'Nf1', 'PDGFB'))                   

                       
pdf(paste0(projDir,'Plots/Cellphonedb.dotplot.wnt.lr.sig.min.2.in.each.row.pdf'), 
        height = 4, 
        width = 6)
                       
## Generate ggplot dotplot
p <- ggplot(data = comb.df4, aes(x = cell.pair, y = interacting_pair, color = sig.genotype, size = mut_wt_neglogpval_fishertest, fill = sig.genotype)) +
  geom_point(shape = 21, aes(colour = as.factor(signif)), alpha = 0.7) +
  scale_colour_manual(values = c("00FFFFFF", "gray", "black")) +  # Outline color
  scale_fill_manual(values = c(genotype.color, "black")) +  # Fill colors
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

# Print the plot
print(p)

# Close the PDF device
dev.off()
                       
print('Done')



##################################################################################
##################################################################################
##################################################################################
### Figure 1e
celltype.name <- 'Tumor'
metaGroupName4 = 'groupID'
genotypes <- unique(srt@meta.data[,metaGroupName4])
# Create a pairwise comparison of genotypes
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
    #df1 <- get(load(input_file_name))
    pair_compare <- paste0(genotype1,'_vs_',genotype2)

    i <- celltype.name
    df <- get(load(input_file_name))
      df <- lig_rec_differential_highlvl_df_not_allzeros[
        grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
        grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE), 
      ]
    df$signif <- 0
    df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
    df$signif <- factor(df$signif, levels = c(0,1))
    df1 <- df
    select_best_pvalue = df[df$mut_wt_neglogpval_fishertest > 0.99,] 

    high_freq_neglogpval_fishertest = names (table (select_best_pvalue$interacting_pair))
    df = df[df$interacting_pair %in% unique(c(high_freq_neglogpval_fishertest)), ]


    df$genotype1 <- genotype1
    df$genotype2 <- genotype2
    df$genotype.pair <- pair_compare
    df.list[[pair_compare]] <- df

    df1$genotype1 <- genotype1
    df1$genotype2 <- genotype2
    df1$genotype.pair <- pair_compare
    df.list1[[pair_compare]] <- df1

}


# List of data frame names you want to rbind
df_names <- c(names(df.list))
df_names1 <- c(names(df.list1))

# Use Reduce with rbind to combine them
comb.df <- Reduce(function(x, y) rbind(x, df.list[[y]]), df_names, init = NULL)
dim(comb.df)
comb.df1 <- Reduce(function(x, y) rbind(x, df.list1[[y]]), df_names1, init = NULL)
dim(comb.df1)

comb.df4 <- comb.df1

# Convert the 'signif' column to numeric
comb.df$signif <- as.numeric(as.character(comb.df$signif))
## selecte only singnificant row 
selected_rows <- comb.df[comb.df$signif >= 1, ]
comb.df3 <- selected_rows


### add higher genotype name in new column
for(n in 1:dim(comb.df3)[1]){
    if(comb.df3[n, "mut_wt_prop_diff"] < 0){
        comb.df3[n, "sig.genotype"] <- comb.df3[n, "genotype2"]
    } else if (comb.df3[n, "mut_wt_prop_diff"] > 0){
        comb.df3[n, "sig.genotype"] <- comb.df3[n, "genotype1"]
    } else if (comb.df3[n, "mut_wt_prop_diff"] == 0){
        comb.df3[n, "sig.genotype"] <- 'No.prop.diff'
    }
}
unique.gene.list <- unique(comb.df3$interacting_pair)                   


# Create an empty dataframe to store the filtered rows
filtered_df2 <- data.frame()

for (un.gene.list.name in unique.gene.list){
    # Filter the dataframe for the current pair
    subset_df <- comb.df3[comb.df3$interacting_pair == un.gene.list.name, ]
    if (length(unique(subset_df$sig.genotype)) == 1){
      filtered_df2 <- rbind(filtered_df2, subset_df)
    }
}

### Get the length of cellpair and ligand receptor pair
interaction.slct <- unique(filtered_df2$interacting_pair)
length(interaction.slct)

un_cell.pair <- unique(filtered_df2$cell.pair)
length(un_cell.pair)

### add higher genotype name in new column
### filter only cell pair and interacting_pair filtered_df2 dataframe
comb.df4 = comb.df4[comb.df4$cell.pair %in% unique(c(un_cell.pair)), ]
comb.df4 = comb.df4[comb.df4$interacting_pair %in% unique(c(interaction.slct)), ]
### add higher genotype name in new column
for(n in 1:dim(comb.df4)[1]){
    if(comb.df4[n, "mut_wt_prop_diff"] == 0.0){
        comb.df4[n, "sig.genotype"] <- 'No.prop.diff'
    } else if (comb.df4[n, "mut_wt_prop_diff"] > 0.0){
        comb.df4[n, "sig.genotype"] <- comb.df4[n, "genotype1"]
    } else if (comb.df4[n, "mut_wt_prop_diff"] < 0.0){
        comb.df4[n, "sig.genotype"] <- comb.df4[n, "genotype2"]
    }
}

new.list <- unique(comb.df4$sig.genotype)

# Move 'No.prop.diff' to the end of the list
levels_sequence <- c(new.list[new.list != 'No.prop.diff'], 'No.prop.diff')
# Convert sig.genotype to a factor with the desired levels sequence
comb.df4$sig.genotype <- factor(comb.df4$sig.genotype, levels = levels_sequence)

# Set up the PDF dimensions
pdf(paste0(projDir,'Plots/Cellphonedb.dotplot.wnt.lr.sig.min.2.in.each.row.pdf'), 
    height = 4, 
    width = 6)

# Convert 'signif' to numeric
comb.df4$signif <- as.numeric(as.character(comb.df4$signif))

### Color
egfr.color <- "#0091CA"
nf1.color <- "#D8423D"
pdgfb.color <- "#55AB55"
genotype.color <- c(egfr.color, nf1.color, pdgfb.color)
#genotype.color <-paletteer_d("ggsci::default_aaas")
# genotype.color <- paletteer_c("grDevices::rainbow", dim(table(comb.df4$sig.genotype)))

#comb.df4$sig.genotype = factor(comb.df4$sig.genotype, levels =unique(comb.df4$sig.genotype)) # set factor #level
comb.df4$sig.genotype = factor(comb.df4$sig.genotype, levels =c('EGFRvIII', 'Nf1', 'PDGFB')) # set factor level                 
gene.name <- 'WNT'
# Filter rows where cell.a or cell.b contains any of the selected cell types
comb.df4 <- comb.df4[grepl(gene.name, comb.df4$gene_a) | grepl(gene.name, comb.df4$gene_b), ] 

# Subsetting based on a single condition
comb.df4 <- comb.df4[comb.df4$genotype.pair != 'Nf1_vs_EGFRvIII', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Excitatory_Neurons.Tumor', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Interneurons.Tumor', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Microglia.Tumor', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Monocytes.Tumor', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Neutrophils.Tumor', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Tumor.Astrocyte', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Tumor.Interneurons', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Tumor.oligodendrocytes', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'Tumor.pDC', ]
comb.df4 <- comb.df4[comb.df4$cell.pair != 'oligodendrocytes.Tumor', ]
comb.df4 <- comb.df4[comb.df4$interacting_pair != 'FZD9_LRP6_WNT5B', ]
comb.df4 <- comb.df4[comb.df4$interacting_pair != 'FZD9_LRP6_WNT5A', ]
comb.df4 <- comb.df4[comb.df4$interacting_pair != 'FZD9_LRP6_WNT3', ]
comb.df4 <- comb.df4[comb.df4$interacting_pair != 'FZD3_LRP6_WNT3', ]
comb.df4 <- comb.df4[comb.df4$interacting_pair != 'FZD2_LRP6_WNT5B', ]
comb.df4 <- comb.df4[comb.df4$interacting_pair != 'FZD2_LRP6_WNT5A', ]

p <- ggplot(data = comb.df4, aes(x = cell.pair, y = interacting_pair, color = sig.genotype, size = mut_wt_neglogpval_fishertest, fill = sig.genotype)) +
  geom_point(shape = 21, aes(colour = as.factor(signif)), alpha = 0.7) +
  scale_colour_manual(values = c("00FFFFFF", "black")) +  # Outline color
  scale_fill_manual(values = c(genotype.color, "black")) +  # Fill colors
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

# Print the plot
print(p)

# Close the PDF device
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure S3a - Wnt gene expression - Muscat analysis
## Here we we used the ourput of muscat analysis result. We ran the muscat analysis at pairwise comparision.
wnt.gene <- c('Wnt5a', 'Lrp5', 'Ptprk', 'Ror2', 'Lrp6', 'Ryk', 'Fzd6', 'Fzd8', 'Fzd7', 'Fzd4', 'Ror1', 'Frzb', 'Epha7', 'Fzd3', 'Wnt7a', 'Wnt7b')
# Variables to specify
srt = srt
#force = FALSE
#do.fgsea = TRUE
logfcThreshold = .5
pvalAdjTrheshold = 0.05
#ds_method = "DESeq2" #c("edgeR", "DESeq2", "limma-trend", "limma-voom")
metaGroupName1 = 'sampleID'
metaGroupName2 = 'celltype'
metaGroupName4 = 'groupID' # Column name in which we want to do pair wise comparision
metaGroupName5 = 'sub_celltype' ## if not assign NULL OR metaGroupName2
org <- 'mouse'
projDir1 <- paste0(projDir,'Plots/Muscat_analysis/')
dir.create(projDir1)
### celltype name to filter column names
celltype.slct <- c('MDM', 'Microglia', 'Monocytes', 'Tcells', 'NKcells', 'Neutrophils', 'Tumor')
gene.list.pathway = list()
# Define the genotypes for comparison
genotypes <- unique(srt@meta.data[,metaGroupName4])
# Create a pairwise comparison of genotypes
genotype_pairs <- combn(genotypes, 2, simplify = FALSE)

new_pair_list = genotype_pairs[[1]]


for (new_pair_list in genotype_pairs) {
    print(genotype1)
    print(genotype2)
    
    
    n = 0
    for (group_name in new_pair_list){
        if (n==0){
            sub_srt <- srt[,srt@meta.data[metaGroupName4] == group_name]
        } else{
            sub_srt1 <- srt[, srt@meta.data[metaGroupName4] == group_name]
            sub_srt <- merge(sub_srt, sub_srt1)
        }
        n = n+1
    }
    sub_srt$muscat_group = ifelse (sub_srt@meta.data[,metaGroupName4] == genotype1,genotype1,genotype2)
    print(table(sub_srt$muscat_group))

    muscatIdents = c(names(table(sub_srt$muscat_group))[1],names(table(sub_srt$muscat_group))[2])
    topGenes = 20
    
    
    pair_compare <- paste0(genotype1,'_vs_',genotype2)
    output_dir <- paste0('../../Genotype.prj/data/Muscat.file/', pair_compare,'/')
    

    tbl_fil1 = readRDS (paste0(output_dir, 'DEGresults_filtered.rds'))
    
    
    # Filter non-significant genes
      tbl_fil2 = lapply(tbl_fil1, function(u) {
      u = dplyr::filter(u, p_adj.loc < pvalAdjTrheshold, abs(logFC) > logfcThreshold)
      dplyr::arrange(u, p_adj.loc)
    })

    tbl_df = do.call (rbind, tbl_fil1)
    colnames (tbl_df)[colnames (tbl_df) == 'logFC'] = 'avg_log2FC'
    colnames (tbl_df)[colnames (tbl_df) == 'p_adj.loc'] = 'p_val_adj'
    colnames (tbl_df)[colnames (tbl_df) == 'cluster_id'] = 'cluster'

    ### Output directory
    out.fig <- projDir1
    #dir.create(out.fig)

    #folder.name <- paste0(genotype1, '_vs_', genotype2, '_cpdb_ds_analysis/')
    ### selecte genes from cellphondb result
#     addGene <- selected_ligand_recptor_from_cellphonedb(cellphonedb.dir, folder.name, metaGroupName4, celltype.slct, out.fig)
#     saveRDS(addGene, paste0(out.fig,'gene.list.rds'))
#     gene.list.pathway[[pair_compare]] <- paste0(out.fig,'gene.list.rds')
    

    addGene = wnt.gene
    
    ## filter our selected genes
    tbl_df <- tbl_df |> filter(gene %in% addGene)

    # Rearrange matrix to Make the heatmap
    avglog2fc.df <- as.matrix(reshape::cast(tbl_df, 
                                            cluster ~ gene, 
                                            value = "avg_log2FC"))

    pvals.df <- as.matrix(reshape::cast(tbl_df, 
                                        cluster ~ gene, 
                                        value = "p_val"))

    # avglog2fc.df <- as.data.frame(avglog2fc.df)
    pvals.df <- as.data.frame(pvals.df)

    # strict cutoffs
    pvals.df[pvals.df <= 0.01] <- "**"
    pvals.df[pvals.df > 0.01 & pvals.df <= 0.05] <- "*"
    pvals.df[is.na(pvals.df)] <- ""
    pvals.df[pvals.df > 0.05] <- ""

    avglog2fc.df[, names(which(sapply(avglog2fc.df, function(x)all(is.na(x)))))] <- NULL # remove cols with all NAs
    pvals.df <- pvals.df[,colnames(avglog2fc.df)] # remove cols with all NAs
    
                                      
    # Calculate the range of values in the matrix 'mat'
    mat_range <- range(avglog2fc.df, na.rm = TRUE)

    # Generate the breaks using the calculated range
    breaks_seq <- seq(mat_range[1], mat_range[2], length.out = 100)
    p <- pheatmap::pheatmap(mat = avglog2fc.df,
                            cluster_rows = F,
                            cluster_cols = F,
                            display_numbers = pvals.df, 
                            breaks = breaks_seq,
                            cellheight = 20, 
                            cellwidth = 20, 
                            color = colorRampPalette(c("navy", "white", "firebrick"))(100))
    pdf(file = paste0(out.fig, pair_compare,
                     ".Wnt.selected.gene.heatmap.muscat.selectedGenes.pdf"), 
        useDingbats = F, 
        height = 2 + ceiling(nrow(avglog2fc.df))/4, 
        width = 4 + ceiling(ncol(avglog2fc.df))/3.5)
    print(p)
    dev.off()
    
}



##################################################################################
##################################################################################
##################################################################################
### Figure 1f and S3b - Wnt gene expression
wnt.gene <- c('Wnt5a', 'Wnt7b', 'Wnt7a', 'Ror1', 'Ror2', 'Lrp5', 'Lrp6', 'Ryk', 'Frzb', 'Ptprk', 'Epha7', 'Fzd3', 'Fzd4', 'Fzd6', 'Fzd7', 'Fzd8')
srt1 <- srt

for (modules_name in wnt.gene){
    #print(modules_name)
    #gene_list <- subset(pathways, term == modules_name)
    modulesL = list (modules_name)
    #print(modulesL)
    message ('Run AddModuleScore')
    srt1 = AddModuleScore (srt1, modulesL, name = modules_name)
}

# Define old and new column names
old_names <- paste0(wnt.gene, '1')
new_names <- wnt.gene
# Get the current column names
current_names <- names(srt1@meta.data)
# Replace old names with new names
current_names[current_names %in% old_names] <- new_names
# Assign the updated column names back to the data frame
names(srt1@meta.data) <- current_names


meta_modules_names = wnt.gene

metaGroupNames = c('orig.ident','groupID')

### Generate boxplots per meta groups
ccomp_df = srt1@meta.data[,meta_modules_names, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt1@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar
color_bar <- color.list[['groupID']]

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
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold"))+
    theme(plot.title = element_text(size = 8)))
      

pdf (paste0 (projDir,'Plots/Wnt.lr_all_genes_module_scores_genes_at_groupID_only.boxplot.seperated.expression.pdf'), width = 14, height = 6)
print (wrap_plots (box_p, ncol= ifelse (length(box_p) > 8,ceiling(length(box_p)/2),length(box_p))))
dev.off()


##################################################################################
##################################################################################
##################################################################################
### Figure 1g and S3c - Wnt gene expression  -TCGA data
wnt.human.gene <- c('WNT5A', 'WNT7B', 'WNT7A', 'ROR1', 'ROR2', 'LRP5', 'LRP6', 'RYK', 'FRZB', 'PTPRK', 'EPHA7', 'FZD3', 'FZD4', 'FZD6', 'FZD7', 'FZD8')
# Create an empty list
gene_list <- list()

gene_list[['WNT.genes.Ligand.recptor']] <- wnt.human.gene
for(gene.name in wnt.human.gene){
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
file.name <- paste0(projDir,"Plots/TCGA.plot_WNT.genes.Ligand.recptor.at_NF1.EGFR.4q12_PDGFRA.pdf")
width = 10
height = 12
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)


##################################################################################
##################################################################################
##################################################################################
### Figure 1h and k - Neftel program expressionat genotype level, presented in quadrant scatter plot
path <- paste0(projDir,'Plots/Neftal_state_plot/')
dir.create(path)
Nefstates1 <- Nefstates(seurat_obj = srt, mouse=T)
color.sch <- paletteer_dynamic("cartography::multi.pal", 20)[c(1, 5, 4, 2)]

p2 <- list()
for (gn in unique(Nefstates1$groupID))
{
    srts1 = subset (Nefstates1, groupID == gn)
    
    p2[[gn]] = ggplot(srts1, aes(x=x_axis, y=y_axis, fill=Tstate, color=Tstate)) +
      geom_point(alpha=0.8, shape=21, size=1.5, color = "black") +
      theme_minimal() +
      xlim (c(min(srts1$x_axis)-0.2, max(srts1$x_axis)+0.2)) +
      ylim (c(min(srts1$y_axis)-0.2, max(srts1$y_axis)+0.2)) +
      xlab ('Relative module score \n[log2(|SC1-SC2|+1)]') +
      ylab ('Relative module score \n[log2(|SC1-SC2|+1)]') +#+ scale_fill_viridis() + scale_color_viridis()
      geom_vline(xintercept = 0,linetype = 'dashed') +
      geom_hline(yintercept = 0,linetype = 'dashed') +
      annotate("text", x = -0.7, y = max(srts1$y_axis[srts1$Tstate == 'OPC-like'])+0.2, label = paste0 (round (sum(srts1$Tstate == "OPC-like")/length(srts1$Tstate) * 100, digits = 2),'%')) +
      annotate("text", x = 0.8, y = max(srts1$y_axis[srts1$Tstate == 'NPC-like'])+0.2, label = paste0 (round (sum(srts1$Tstate == "NPC-like")/length(srts1$Tstate) * 100, digits = 2),'%')) +
      annotate("text", x = -0.7, y = min(srts1$y_axis[srts1$Tstate == 'AC-like'])-0.2, label = paste0 (round (sum(srts1$Tstate == "AC-like")/length(srts1$Tstate) * 100, digits = 2),'%')) +
      annotate("text", x = 0.8, y = min(srts1$y_axis[srts1$Tstate == 'MES-like'])-0.2, label = paste0 (round (sum(srts1$Tstate == "MES-like")/length(srts1$Tstate) * 100, digits = 2),'%')) +
      ggtitle(gn)+
      #scale_color_manual(values=c("black", "red", "green","blue")) +
    scale_color_manual(values=color.sch, guide = guide_legend(override.aes = list(size = 3))) +
      scale_fill_manual(values=color.sch, guide = guide_legend(override.aes = list(size = 3)))+ theme(axis.title=element_text(size=13)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", color = "black"),  # Set x-axis text to bold
        axis.text.y = element_text(face = "bold", color = "black"),  # Set y-axis text to bold
        axis.title = element_text(face = "bold", color = "black"),  # Set axis titles to bold
        strip.text = element_text(face = "bold", color = "black"),  # Set facet strip text to bold
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
         )

    png (paste0(path, 'Diagram_GBM_states_',gn,'_Nefstates_function1.png'), height=1100, width=1400, res=300)
    print (p2[[gn]])
    dev.off()

}


### Cell cycle plot

p2 <- list()
p3 <- list()
for (gn in unique(Nefstates1$groupID))
{
    srts1 = subset (Nefstates1, groupID == gn)

    filtered_data <- srts1[srts1$Tstate == "OPC-like", ]
    OPC.cc = round (sum(filtered_data$Phase2 == "cycling")/length(filtered_data$Phase2) * 100, digits = 2)
    OPC.non.cc = round (sum(filtered_data$Phase2 == "non-cycle")/length(filtered_data$Phase2) * 100, digits = 2)

    filtered_data <- srts1[srts1$Tstate == "NPC-like", ]
    NPC.cc = round (sum(filtered_data$Phase2 == "cycling")/length(filtered_data$Phase2) * 100, digits = 2)
    NPC.non.cc = round (sum(filtered_data$Phase2 == "non-cycle")/length(filtered_data$Phase2) * 100, digits = 2)

    filtered_data <- srts1[srts1$Tstate == "AC-like", ]
    AC.cc = round (sum(filtered_data$Phase2 == "cycling")/length(filtered_data$Phase2) * 100, digits = 2)
    AC.non.cc = round (sum(filtered_data$Phase2 == "non-cycle")/length(filtered_data$Phase2) * 100, digits = 2)

    filtered_data <- srts1[srts1$Tstate == "MES-like", ]
    MES.cc = round (sum(filtered_data$Phase2 == "cycling")/length(filtered_data$Phase2) * 100, digits = 2)
    MES.non.cc = round (sum(filtered_data$Phase2 == "non-cycle")/length(filtered_data$Phase2) * 100, digits = 2)
    
    
#     cycling_color <- "#ff7f0e"  # Orange for cycling cells
#     non_cycling_color <- "#CCCCCC"  # Light green for non-cycling cells
cycling_color <- "#006400"  # Dark green for cycling cells
non_cycling_color <- "#ffc0cb"  # Light pink for non-cycling cells

    cc.color <- c(cycling_color, non_cycling_color)
    
    
    p3[[gn]] = ggplot(srts1, aes(x=x_axis, y=y_axis, fill=cc_score, color=cc_score)) +
    #geom_point(alpha=0.8, shape=21, size=0.1) +
    geom_point(alpha=0.8, shape=21, size=0.5) +
    #scale_fill_viridis(option='rocket') + scale_color_viridis(option='rocket') +
    scale_fill_gradient(low = non_cycling_color, high = cycling_color) +  scale_color_gradient(low = non_cycling_color, high = cycling_color) +
    theme_minimal() +
    xlim (c(min(srts1$x_axis)-0.2, max(srts1$x_axis)+0.2)) +
    ylim (c(min(srts1$y_axis)-0.2, max(srts1$y_axis)+0.2)) +
    xlab ('Relative module score \n[log2(|SC1-SC2|+1)]') +
    ylab ('Relative module score \n[log2(|SC1-SC2|+1)]') +#+ scale_fill_viridis() + scale_color_viridis()
    geom_vline(xintercept = 0,linetype = 'dashed') +
    geom_hline(yintercept = 0,linetype = 'dashed') +
    annotate("text", x = -0.7, y = max(srts1$y_axis[srts1$Tstate == 'OPC-like'])+0.2, label = paste0 ("OPC - CC", '(',OPC.cc,'%)'), size = 3) +
    annotate("text", x = 0.8, y = max(srts1$y_axis[srts1$Tstate == 'NPC-like'])+0.2, label = paste0 ("NPC - CC", '(',NPC.cc,'%)'), size = 3) +
    annotate("text", x = -0.7, y = min(srts1$y_axis[srts1$Tstate == 'AC-like'])-0.2, label = paste0 ("AC - CC", '(',AC.cc,'%)'), size = 3) +
    annotate("text", x = 0.8, y = min(srts1$y_axis[srts1$Tstate == 'MES-like'])-0.2, label = paste0 ("MES - CC", '(',MES.cc,'%)'), size = 3) +
    ggtitle(gn) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),  # Set x-axis text to bold
        axis.text.y = element_text(face = "bold"),  # Set y-axis text to bold
        axis.title = element_text(face = "bold"),  # Set axis titles to bold
        strip.text = element_text(face = "bold"),  # Set facet strip text to bold
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.line = element_line(color = "black")
         )
    
    png (paste0(path, 'Diagram_GBM_states_',gn,'_cc_score.png'), height=1100, width=1400, res=300)
    print (p3[[gn]])
    dev.off()
}


### Figure 1i - Neftel program expression genotype quantify and presented in boxplot
Nefstates1$Tstate1[Nefstates1$Tstate %in% c('MES-like', 'AC-like')] = 'MES.AC-like'
Nefstates1$Tstate1[Nefstates1$Tstate %in% c('NPC-like', 'OPC-like')] = 'NPC.OPC-like'
srt$Tstate1 <- Nefstates1$Tstate1[match(colnames(srt), rownames(Nefstates1))]

## Parameter
metaGroupName1 = 'orig.ident'
metaGroupName2 = 'groupID'
metaGroupName3 = 'Tstate1'
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


pdf (paste0 (path, 'cell_composition_',metaGroupName3,'_Group_id_level_barboxplots_t_test_with_padj.pdf'), width=8, height=5)
(cc_bar | p) + plot_layout (widths= c(1,celltype_length))
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure 1j - Neftel program expression in human data - quantify and presented in boxplot
library ('clusterProfiler')
gbm_modules_6states = readRDS (paste0 ('../../Genotype.prj/data/human_GBM_netfel_6_states.rds'))
gene_list <- list()
#gbm_modules_6states
gene_list[['MES.AC']] <- c(gbm_modules_6states['MES1'][[1]], gbm_modules_6states['MES2'][[1]], gbm_modules_6states['AC'][[1]])
gene_list[['NPC.OPC']] <- c(gbm_modules_6states['NPC1'][[1]], gbm_modules_6states['NPC2'][[1]], gbm_modules_6states['OPC'][[1]])

result <- Add_celltype_expression(gbm.2018, gene_list)
tmp <- result[[1]]
modules <- result[[2]]

EGFR.color <- "#0091CA"
NF1.color <- "#D8423D"
PDGFB.color <- "#55AB55"
genotype.color <- c(EGFR.color, NF1.color, PDGFB.color)
col.name <- 'new_subtype_based_on_mut_selected_from_supp_7_subtype'
level.factor <- c("EGFRvIII", "NF1", "4q12_PDGFRA")
file.name <- paste0(path,"TCGA_gbm.Neftel.prog_at_NF1.EGFR.PDGFB.pdf")
width = 5
height = 5
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)


##################################################################################
##################################################################################
##################################################################################
### Figure S4 - Neftel prog UMAP
gbm_modules_6states = readRDS (paste0 ('../../Genotype.prj/data/mouse_GBM_netfel_6_states.rds'))
gbm_modules_6states[['MES1']] <- NULL
gbm_modules_6states[['NPC2']] <- NULL
srt1 <- srt
for (modules_name in names(gbm_modules_6states)){
    message (modules_name)
    message ('Run AddModuleScore')
    srt1 = AddModuleScore (srt1, list(gbm_modules_6states[[modules_name]]), name = modules_name)
}
meta_modules_names = paste0(names(gbm_modules_6states), '1') # remove grey module (unassigned)
umap_df = data.frame (srt1[[reductionName]]@cell.embeddings, srt1@meta.data[,meta_modules_names])
umap_p1 = lapply (meta_modules_names, function(x) ggplot(data = umap_df) + 
geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = 0.1) + 
scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
theme_bw() + ggtitle (x))

png (paste0 (path,'Netfel_all_6_prog_module_scores_umap.png'), width = 2500, height = 2500, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 4,ceiling(length(umap_p1)/4),length(umap_p1))))
dev.off()


##################################################################################
##################################################################################
##################################################################################
### Figure 1l - In this figure we are showing the RNA velocity which we calculate the by standerd RNA velocity pipeline by scVelo and we added Tumor compartment embadding with cell cycle information to show the RNA velocity UMAP.




##################################################################################
##################################################################################
##################################################################################
### Figure 1m - WGCNA identified program - Calculate mean module score and presented in boxplot
wgcna_gmt <- GSA.read.gmt(paste0('../../Genotype.prj/data/Tumor.wgcna_14modules_25genes_genesets.gmt.txt'))
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
new_colnames = c('Organogenesis.Ependymal.Cell', 'TNFA.Sig', 'Interferon', 'Myc.Target.Ox.Phos', 
                 'Kras.Sign1', 'COMPLEMENT.Innate.Immune.Sytm',
                 'EMT', 'Histone.Expressing.Cluster', 'Interferon.Gamma.Allograft', 'G2M', 'Hypoxia.Glycolysis', 'Kras.Sign2', 'G1S', 'Myc.Target')
colnames(srt1@meta.data)[(length(colnames(srt@meta.data))+1):length(colnames(srt1@meta.data))] <- new_colnames
metaGroupNames = c('orig.ident','groupID')
srt_wgcna <- srt1
projDirW <- paste0(projDir, 'Plots/WGCNA.analysis/')
dir.create(projDirW)


### Make plot with selected module
new_colnames1 = c('TNFA.Sig', 'Interferon',  'G2M')
meta_modules_names = new_colnames1
metaGroupNames = c('orig.ident','groupID')
ccomp_df = srt_wgcna@meta.data[,new_colnames1, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar
color_bar <- color.list[[metaGroupNames[2]]]
box_p = lapply (seq_along(meta_modules_names), function(x) 
  ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
    geom_boxplot(fill = color_bar, outlier.shape = NA) +
    geom_jitter (color="black", size=0.4, alpha=0.9) +
    geom_pwc(method = "wilcox_test", label = "{p.format}{p.signif}", hide.ns = T, p.adjust.method = "none") +
    theme_bw() + 
    scale_fill_manual (values= color_bar) + 
    ggtitle (meta_modules_names[x]) + 
    theme_classic()+
               theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", color = "black"),          # Set x-axis text to bold
        axis.text.y = element_text(face = "bold", color = "black"),  # Set y-axis text to bold
        axis.title = element_text(face = "bold", color = "black"),  # Set axis titles to bold
        strip.text = element_text(face = "bold", color = "black"),  # Set facet strip text to bold
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(face = "bold", color = "black", size = 12)
         ))

pdf (paste0 (projDirW,'WGCNA_module_scores_umap_new_ann_at_groupID_only.boxplot.selected.modules.top25.genes.pdf'), width = 6, height = 3)
print (wrap_plots (box_p, ncol= ifelse (length(box_p) > 8,ceiling(length(box_p)/2),length(box_p))))
dev.off()




##################################################################################
##################################################################################
##################################################################################
### Figure 1m - WGCNA identified program - Calculate mean module score and presented in boxplot
load(file = paste0("../../Genotype.prj/data/Human.converted.wgcna_14modules_25genes_genesets.gmt.Rda"))
result <- Add_celltype_expression(gbm.2018, gene_list.prog1[c('Interferon', 'TNFA.Sig', 'G2M')])
tmp <- result[[1]]
modules <- result[[2]]
EGFR.color <- "#0091CA"
NF1.color <- "#D8423D"
PDGFB.color <- "#55AB55"
genotype.color <- c(EGFR.color, NF1.color, PDGFB.color)
col.name <- 'new_subtype_based_on_mut_selected_from_supp_7_subtype'
level.factor <- c("EGFRvIII", "NF1", "4q12_PDGFRA")
file.name <- paste0(projDirW,"TCGA.DATA.WGCNA.PROG.FROM.selected.program_at_NF1.EGFR.4q12_PDGFRA.pdf")
width = 4
height = 6
boxplot.fun1(obj = tmp, col.name = col.name, genotype.color = genotype.color, level.factor = level.factor, file.name = file.name,
                       width = width, height = height)


##################################################################################
##################################################################################
##################################################################################
### Figure S4c  - Hallmark programm expression
library ('clusterProfiler')
gmt.file = paste0 ('../../Genotype.prj/data/gsea.hallmark.prog.mh.all.v2023.1.Mm.symbols.gmt')
pathways = read.gmt(gmt.file)
pathways <- subset(pathways, term == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB' | term == 'HALLMARK_INTERFERON_ALPHA_RESPONSE' | term == 'HALLMARK_INTERFERON_GAMMA_RESPONSE' | term == 'HALLMARK_G2M_CHECKPOINT')

srt1 <- srt
for (modules_name in unique(pathways$term)){
    #print(modules_name)
    gene_list <- subset(pathways, term == modules_name)
    modulesL = list (gene_list$gene)
    #print(modulesL)
    message ('Run AddModuleScore')
    srt1 = AddModuleScore (srt1, modulesL, name = modules_name)
}

metaGroupNames = c('orig.ident','groupID')
meta_modules_names = paste0(unique(pathways$term), '1')
### Generate boxplots per meta groups
ccomp_df = srt1@meta.data[,meta_modules_names, drop=FALSE]
#ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 
ccomp_df = aggregate (ccomp_df, by=as.list(srt1@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar
color_bar <- color.list[[metaGroupNames[2]]]

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
    #facet_wrap (as.formula(paste("~", metaGroupNames[2]))) 
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"))+
    theme(plot.title = element_text(size = 8)))
#}        

pdf (paste0 (projDir,'Plots/Hallmark_gene_module_scores_genes_at_groupID_only.boxplot.only.emt.inter.g2m.tnfa.pdf'), width = 12, height = 4)
print (wrap_plots (box_p, ncol= ifelse (length(box_p) > 7,ceiling(length(box_p)/8),length(box_p))))
dev.off()



##################################################################################
##################################################################################
##################################################################################
### Figure S8b - EGFR expression in mouse scRNA data
srt1 <- srt
modulesL = list ('Egfr')
message ('Run AddModuleScore')
srt1 = AddModuleScore (srt1, modulesL, name = 'Egfr')

new_colnames1 <- c('Egfr1')
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

pdf (paste0 (projDir,'Plots/Egfr_module_scores.boxplot.pdf'), width = 3, height = 3)
print (wrap_plots (box_p, nrow= 1))
dev.off()


