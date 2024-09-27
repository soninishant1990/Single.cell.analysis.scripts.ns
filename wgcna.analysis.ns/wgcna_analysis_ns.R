


# load libraries
library('dplyr')
library('Seurat')
library('patchwork')
library('ggplot2')
library('harmony')
library('paletteer')
library ('RColorBrewer')
library('hash')
library ('WGCNA')
library ('hdWGCNA')
library ('igraph')
library ('viridis')
library ('clusterProfiler')
library('tidyr')
library('ggpubr')

# ### Load WGCNA source scripts
# wgcna_path <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_3_pediatric/LabCode/hdWGCNA-dev/R/'
# script_files <- list.files(path = wgcna_path, pattern = "*.R", full.names = TRUE)
# for (script_file in script_files) {
#   source(script_file)
# }

allowWGCNAThreads()
projDirW = paste0(wgcna_dir,'WGCNA_sp_',softPower,'_ds_',deepSplit,'_mch_',mergeCutHeight,'max_shared_',max_shared,'_minModuleSize_',minModuleSize,'_metasize_',metacells_k,'_genes_',length(genes.keep),'/')
system (paste('mkdir -p',paste0(projDirW,'Plots/')))
setwd (projDirW)

if (!file.exists (paste0(projDirW, 'srt_wgcna.rds')) | force) {
    message ('No suerat object found!')
    srt_wgcna <- SetupForWGCNA(
    srt_wgcna,
    gene_select = "custom", # the gene selection approach
    features = genes.keep,  
    #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "hdWGCNA" # the name of the hdWGCNA experiment
    )


    # construct metacells  in each group
    filter_sample = table (srt_wgcna@meta.data[,metacells_groups[1]])
    filter_sample = names (filter_sample)[filter_sample > metacells_k]
    filter_cells = colnames(srt_wgcna) [srt_wgcna@meta.data[,metacells_groups[1]] %in% filter_sample]

    srt_wgcna = srt_wgcna[, colnames (srt_wgcna) %in% filter_cells]
    srt_wgcna <- MetacellsByGroups (
    seurat_obj = srt_wgcna,
    group.by = unique(metacells_groups), # specify the columns in seurat_obj@meta.data to group by
    k = metacells_k, # nearest-neighbors parameter
    max_shared = max_shared,
    ident.group = metacells_groups[1] # set the Idents of the metacell seurat object
    )

    # Number of metacells
    message (paste ('number of metacells:', ncol(GetMetacellObject(srt_wgcna))))
    # normalize metacell expression matrix:
    srt_wgcna <- NormalizeMetacells (srt_wgcna)

    srt_wgcna <- SetDatExpr(
    srt_wgcna,
    group_name = unique(srt_wgcna@meta.data[,metacells_groups[2]]), # the name of the group of interest in the group.by column
    group.by=metacells_groups[2] # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    )


    if (powerTable)
    {
    srt_wgcna <- TestSoftPowers (
    srt_wgcna,
    powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
    setDatExpr = FALSE # set this to FALSE since we did this above
    )
    plot_list <- PlotSoftPowers(srt_wgcna)

    pdf(paste0(projDirW, "Plots/Power_threshold_plots.pdf"), height=7, width=8)
    print (wrap_plots (plot_list))
    dev.off()
    }


    # construct co-expression network:
    srt_wgcna <- ConstructNetwork (
    srt_wgcna, 
    soft_power=softPower,
    setDatExpr=FALSE,
    min_power = 3,
    tom_outdir = paste0(projDirW,"TOM"),
    use_metacells = TRUE,
    group.by = NULL,
    group_name = NULL,
    consensus = FALSE,
    multi.group.by = NULL,
    multi_groups = NULL,
    blocks = NULL,
    maxBlockSize = 30000,
    randomSeed = 12345,
    corType = "pearson",
    consensusQuantile = 0.3,
    networkType = "signed",
    TOMType = "signed",
    TOMDenom = "min",
    scaleTOMs = TRUE,
    scaleQuantile = 0.8,
    sampleForScaling = TRUE,
    sampleForScalingFactor = 1000,
    useDiskCache = TRUE,
    chunkSize = NULL,
    deepSplit = deepSplit,
    pamStage = FALSE,
    detectCutHeight = 0.995,
    minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    saveConsensusTOMs = TRUE,
    #tom_name = 'all',
    overwrite_tom=TRUE
    )
    #srt_wgcna@misc$hdWGCNA$wgcna_net$TOMFiles = paste0(projDirW,'TOM/_TOM.rda')
    srt_wgcna@misc$hdWGCNA$wgcna_net$TOMFiles = paste0(projDirW,'TOM/hdWGCNA_TOM.rda')

    pdf(paste0(projDirW,'Plots/WGCNA_SignedDendro.pdf'),height=5, width=8)
    print (PlotDendrogram(srt_wgcna, main='hdWGCNA Dendrogram'))
    dev.off()

    # need to run ScaleData first or else harmony throws an error:
    srt_wgcna <- ScaleData (srt_wgcna, features=VariableFeatures(srt_wgcna))


    # compute all MEs in the full single-cell dataset
    srt_wgcna <- ModuleEigengenes(
    srt_wgcna,
    group.by.vars=metaGroupNames[2]
    )
    modules <- GetModules(srt_wgcna)

    #ComputeModuleEigengene
    srt_wgcna <- ModuleConnectivity(
    srt_wgcna,
    group.by = NULL, 
    group_name = NULL,
    wgcna_name = NULL
    )

    #srt_wgcna@misc$hdWGCNA$wgcna_net$TOMFiles = paste0(projDirW,'TOM/hdWGCNA_TOM.rda')

    srt_wgcna <- RunModuleUMAP(
    srt_wgcna,
    n_hubs = 10, # number of hub genes to include for the UMAP embedding
    n_neighbors=15, # neighbors parameter for UMAP
    min_dist=0.1 # min distance between points in UMAP space
    )

    ### AddModuleScore of WGCNA modules ###
    modules <- GetModules(srt_wgcna)
    modulesL = split (modules$gene_name, modules$module)
    message ('Run AddModuleScore for each WGCNA module')
    srt_wgcna = ModScoreCor (
        seurat_obj = srt_wgcna, 
        geneset_list = modulesL, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'wgcna_modules', outdir = paste0(projDirW,'Plots/'))

    # Assign cells to modules based on highest modules score
    message ('Hard-assign cells to wgcna modules')
    modulesL = split (modules$gene_name, modules$module)

    # Reduce size of seurat object before saving
    srt_wgcna = DietSeurat(
    srt_wgcna,
    counts = FALSE,
    data = TRUE,
    scale.data = FALSE,
    features = NULL,
    assays = NULL,
    dimreducs = c(reductionName, reductionSave),
    graphs = NULL
    )
    print(table(modules$module))
    message('Saving WGCNA seurat object')
    saveRDS (srt_wgcna, paste0(projDirW, 'srt_wgcna.rds'))

    message('Saving module dataframe with kME value')
    # Assuming you have already loaded the dataframe into a variable called 'df'

    # Selecting 'gene_name' and 'module' columns
    new_df <- modules[, c('gene_name', 'module')]

    # Creating an empty column 'kME_value'
    new_df$kME_value <- NA

    # Looping through each row in the dataframe
    for (i in 1:nrow(new_df)) {
    # Getting the module value for the current row
    module <- new_df[i, 'module']
    
    # Extracting the corresponding kME column name
    kME_col <- paste0('kME_', module)
    
    # Checking if the kME column exists and is not missing
    if (kME_col %in% colnames(modules) && !is.na(modules[i, kME_col])) {
        # Getting the kME value from the matching column
        kME_value <- modules[i, kME_col]
        
        # Assigning the kME value to the 'kME_value' column
        new_df[i, 'kME_value'] <- kME_value
    }
    }


    # Printing the resulting dataframe
    head(new_df)
    ## Save the dataframe
    write.csv (new_df, paste0(projDirW, 'module_assignments.csv'))
    write.csv (print(table(modules$module)), paste0(projDirW, 'module_names.csv'))

} else {
    srt_wgcna <- readRDS(paste0(projDirW, 'srt_wgcna.rds'))
    modules <- GetModules(srt_wgcna)
    modulesL = split (modules$gene_name, modules$module)
}


# Run GSEA enrichment on each cluster DE markers
if (!file.exists (paste0(projDirW, 'EnrichR_WGCNA_module_genes.rds')) | force)
   {
    gmt_annotations = c(
    'h.all.v7.1.symbol.gmt',
    'c2.cp.kegg.v7.1.symbol.gmt',
    'c2.cp.reactome.v7.1.symbol.gmt',#,
    'c5.bp.v7.1.symbol.gmt'#, # GO_BP
    #'c3.tft.v7.1.symbol.gmt' # TF database
    )
    
    # GSEA analysis on DEG per cluster
    EnrichRResAll = list()
    for (ann in gmt_annotations)
      {
      gmt.file = paste0 ('/ahg/regevdata/projects/ICA_Lung/Bruno/DBs/GSEA_gs/',org,'/',ann)
      pathways = read.gmt(gmt.file)
      #pathways = gmtPathways (gmt.file)
      message (paste('Compute enrichment per cluster using annotation:', ann))
      EnrichRResCluster = list()
      for (i in unique (modules$module))
        {
        message (paste ('EnrichR running module',i)) 
        sig_genes = modules$gene_name[modules$module == i]
        egmt <- enricher (sig_genes, TERM2GENE=pathways, universe = enricher_universe)
          if (length(egmt) == 0){
                message (paste ('Zoro genes match with Enricher')) 
            } else{
              EnrichRResCluster[[i]] = egmt@result
          }
        
        }
      EnrichRResAll[[ann]] = EnrichRResCluster
      }
    saveRDS (EnrichRResAll, paste0(projDirW, 'EnrichR_WGCNA_module_genes.rds'))  
   #} else {
    # EnrichRResAll <- readRDS(paste0(projDirW, 'EnrichR_WGCNA_module_genes.rds'))
   #}
    pvalAdjTrheshold = 0.05
    top_pathways = 10
    EnrichRes_dp = lapply (EnrichRResAll, function(x) dotGSEA (enrichmentsTest_list = x, type = 'enrich', padj_threshold = pvalAdjTrheshold, top_pathways= top_pathways))
    
    for (i in seq_along(gmt_annotations))
      {
        if (i == 1){
          width = 7
          hight = 7
        } else if (i==2){
            width = 8
            hight = 8
        } else if (i==3){
            width = 13
            hight = 13
        } else if (i==4){
            width = 8
            hight = 12
        }

      pdf (paste0(projDirW,'Plots/EnrichR_WGCNA_module_genes_',gmt_annotations[[i]],'_dotplots.pdf'),width, hight)
      print (EnrichRes_dp[[i]])
      dev.off()
      
      }  
   }


if (do.plots)
{

    # hubgene network
    message ('Generate gene network plot')
    pdf (paste0(projDirW, 'Plots/hubGeneNetwork.pdf'))
    HubGeneNetworkPlot(
    srt_wgcna,
    n_hubs = 8, n_other=10,
    edge_prop = 0.75,
    mods = 'all'
    )
    dev.off()

    pdf (paste0(projDirW, 'Plots/hubGeneNetwork_UMAP.pdf'))
    ModuleUMAPPlot(
    srt_wgcna,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1, # proportion of edges to sample (20% here)
    label_hubs=2 ,# how many hub genes to plot per module?
    keep_grey_edges=FALSE
    )
    dev.off()

    # modulescore of modules on UMAPs
    meta_modules_names = names(modulesL)[names(modulesL) != 'grey'] # remove grey module (unassigned)
    umap_df = data.frame (srt_wgcna[[reductionName]]@cell.embeddings, srt_wgcna@meta.data[,meta_modules_names])

    umap_p1 = lapply (meta_modules_names, function(x) ggplot(data = umap_df) + 
    geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    theme_bw() + ggtitle (x))

    ### Generate boxplots per meta groups
    ccomp_df = srt_wgcna@meta.data[,names (modulesL), drop=FALSE]
    #ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 
    ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)

    ## Get the color bar
    color_bar <- paletteer_c("grDevices::rainbow", dim(table(srt_wgcna@meta.data[,metaGroupNames[2]])))


    #if (exists ('order_factor')) ccomp_df$Genotype = factor (ccomp_df$Genotype, levels = order_factor)
    if (length(unique (srt_wgcna@meta.data[,metaGroupNames[3]])) == length (unique(srt_wgcna@meta.data[,metaGroupNames[2]])))
    {
    box_p = lapply (seq_along(meta_modules_names), function(x) 
      ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= meta_modules_names[x],fill = metaGroupNames[3])) +
            #geom_violin (trim=TRUE) +
            geom_bar (stat="identity") +
            theme_bw() + 
            scale_fill_manual (values= module_pal) + 
            ggtitle (meta_modules_names[x]) + 
            facet_wrap (as.formula(paste("~", metaGroupNames[1]))) + theme_classic()) 
    } else {
    box_p = lapply (seq_along(meta_modules_names), function(x) 
      ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
        #geom_violin (trim=TRUE) +
        geom_boxplot () +
        geom_boxplot(fill = color_bar) +
        geom_jitter (color="black", size=0.4, alpha=0.9) +
        geom_pwc(method = "wilcox_test", label = "{p.format}{p.signif}", hide.ns = T, p.adjust.method = "none") +
        theme_bw() + 
        scale_fill_manual (values= color_bar) + 
        ggtitle (meta_modules_names[x]) + 
        #facet_wrap (as.formula(paste("~", metaGroupNames[2]))) 
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    }        

    png (paste0 (projDirW,'Plots/WGCNA_module_scores_umap.png'), width = 8000, height = 4000, pointsize=10, res = 300, type="cairo")
    print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))) / 
    wrap_plots (box_p, ncol= ifelse (length(box_p) > 8,ceiling(length(box_p)/2),length(box_p))) + plot_layout(heights=c(1,1.5)))
    dev.off()



    box_p = lapply (seq_along(meta_modules_names), function(x) 
    ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= meta_modules_names[x])) +
      #geom_violin (trim=TRUE) +
      geom_boxplot () +
      geom_boxplot(fill = color_bar) +
      geom_jitter (color="black", size=0.4, alpha=0.9) +
      geom_pwc(method = "wilcox_test", label = "{p.format}{p.signif}", hide.ns = T, p.adjust.method = "none") +
      theme_bw() + 
      scale_fill_manual (values= color_bar) + 
      ggtitle (meta_modules_names[x]) + 
      #facet_wrap (as.formula(paste("~", metaGroupNames[2]))) 
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  #}        

  png (paste0 (projDirW,'Plots/WGCNA_module_scores_umap1.png'), width = 8000, height = 4000, pointsize=10, res = 300, type="cairo")
  print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))) / 
  wrap_plots (box_p, ncol= ifelse (length(box_p) > 8,ceiling(length(box_p)/2),length(box_p))) + plot_layout(heights=c(1,1.5)))
  dev.off()

    # Plot top 25 genes per module
    for (i in unique(modules$module)) {
    message('Generate module networks for top 25 genes')
    if (i != 'gray') {
        tryCatch({
            ModuleNetworkPlot(srt_wgcna, mods = i, outdir = paste0(projDirW, 'ModuleNetworks'))
        }, error = function(e) {
            if (grepl("subscript out of bounds", e$message)) {
                message("Error: subscript out of bounds for module ", i)
            } else {
                message("An unexpected error occurred for module ", i, ": ", e$message)
            }
        })
    }
  }
}
