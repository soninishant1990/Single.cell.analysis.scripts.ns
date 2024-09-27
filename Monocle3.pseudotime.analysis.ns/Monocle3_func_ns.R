

### parameter
## TUMOR
# projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Tumor.comp/Tumor.comp_subset/sampleID_harmony/Tumor.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/Tumor.after.removing.doublet_subset/sampleID_harmony/Tumor.after.regress.out.nfeature_subset/sampleID_harmony/'
# srt <- readRDS(paste0(projDir,'srt.rds'))
# ## Provide the column name for plot learn graph
# metaGroupNames=c('groupID', 'sampleID', 'celltype', 'celltype')
# reductionName = "sampleID_harmony_umap"
# Monocle.dir.name = 'Monocle.with.harmony'
# Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name)



Monocle.fun <- function(srtS = NULL, metaGroupNames = NULL, reductionName = NULL, projDir = NULL, Monocle.dir.name = NULL, org = NULL) {

    message('Create dir for Monocle outputs')
    projDirS <- paste0(projDir, Monocle.dir.name, '/')
    dir.create(paste0(projDirS, 'Plots/'), recursive = TRUE, showWarnings = FALSE)

    # Check if any value from metaGroupNames is not present in colnames(srtS@meta.data)
    if (any(!metaGroupNames %in% colnames(srtS@meta.data))) {
        stop("Code stopped: Some values from metaGroupNames are not present in colnames(srtS@meta.data).")
    }

    message ('Cell cycle analysis')
    if (org == 'Mouse') {
        m.g2m.genes = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/cellcycle_mouse.g2m.genes.rds')
        m.s.genes = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/cellcycle_mouse.s.genes.rds')
        s.genes <- m.s.genes
        g2m.genes <- m.g2m.genes
    } else if (org == 'Human'){
        s.genes <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes
    }
    
    ## Calculate CellCycleScoring score
    srtS <- CellCycleScoring(srtS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

    metaGroupNames1 = c('Phase')
    umap <- DimPlot (object = srtS, reduction = reductionName, pt.size = 1, cols =paletteer_c("grDevices::rainbow", dim(table(srtS@meta.data[,metaGroupNames1]))), group.by = metaGroupNames1, label = F) +theme(legend.position="bottom")
    png (paste0(projDirS, 'Plots/Cell_cycle_phase_umap.png'), width = 2200, height = 2200, pointsize=10, res = 300, type="cairo")
    print (wrap_plots (umap))
    dev.off()

    message('convert seurat to cds obj')
    cds_file <- paste0(projDirS, 'cds.rds')

    if (file.exists(cds_file)) {
    message('Reading existing cds object...')
    cds <- readRDS(cds_file)
    } else {
    # If the file doesn't exist, perform the conversion steps
    cds <- as.cell_data_set(srtS)
    
    ## Calculate size factors using built-in function in monocle3
    cds <- estimate_size_factors(cds)
    
    ## Add gene names into CDS
    cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srtS[["RNA"]])
    cds <- cluster_cells(cds)
    
    ######################################################
    recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
    names(recreate.partitions) <- cds@colData@rownames
    recreate.partitions <- as.factor(recreate.partitions)
    cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
    cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- srtS@reductions[[reductionName]]@cell.embeddings
    
    message('learn graph')
    cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)
    
    message('Plot learn_graph')
    genes <- NULL
    for (i in 1:length(metaGroupNames)) {
        if (length(metaGroupNames[i]) != 1) {
        pdf(paste0(projDirS, 'Plots/Monocle_srt_learn_graph_', dim(cds)[2], '_cells_umaps_', metaGroupNames[i], '.pdf'), 8, 5)
        print(plot_cells(cds,
                        color_cells_by = metaGroupNames[i],
                        label_principal_points = TRUE,
                        label_cell_groups = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = TRUE,
                        graph_label_size = 1.5,
                        label_roots = TRUE))
        if (!is.null(genes)) print(plot_cells(cds,
                                                label_cell_groups = FALSE,
                                                genes = genes,
                                                label_leaves = TRUE,
                                                label_branch_points = TRUE,
                                                graph_label_size = 1.5))
        dev.off()
        }
    }
    
    message('Plot show on show_trajectory_graph')
    p1 <- plot_cells(cds, show_trajectory_graph = TRUE, label_cell_groups = FALSE)
    pdf(paste0(projDirS, 'Plots/Monocle_srt_umaps_harmony.pdf'), width = 8.5, height = 7)
    print(wrap_plots(p1))
    dev.off()
    
    message('Save cds object...')
    saveRDS(cds, cds_file)
    }
  
    for (i in 1:length(metaGroupNames)) {
        metaGroupName <- metaGroupNames[[i]]
        print(metaGroupName)
        pdf(paste0(projDirS, 'Plots/UMAPs', '_', metaGroupName, '_wHarmony_new_rigor.pdf'), width = 8, height = 5)
        umap <- DimPlot(srtS, reduction = reductionName, group.by = metaGroupName, label = FALSE, pt.size = .1, label.size = 5)
        plot(umap)
        dev.off()
    }
  
    message('Plot show Color cells by pseudotime')
    cluster_levels <- as.numeric(cds@clusters$UMAP$clusters)
    max_cluster_level <- max(cluster_levels)
    root5 <- order_cells(cds, root_cells = colnames(cds[, clusters(cds) == max_cluster_level]))
    
    p1 <- plot_cells(root5,
                    color_cells_by = "pseudotime",
                    group_cells_by = "cluster",
                    label_cell_groups = FALSE,
                    label_groups_by_cluster = FALSE,
                    label_leaves = FALSE,
                    label_branch_points = FALSE,
                    label_roots = FALSE,
                    trajectory_graph_color = "black")
    
    pdf(paste0(projDirS, 'Plots/Color_cells_by_pseudotime1.pdf'), width = 7, height = 5)
    print(wrap_plots(p1))
    dev.off()
    
    png(paste0(projDirS, 'Plots/Color_cells_by_pseudotime1.png'), width = 1700, height = 1200, res = 300)
    print(wrap_plots(p1))
    dev.off()
    
    message('Identify genes that change as a function of pseudotime')
    cds_graph_test_results <- graph_test(cds,
                                        neighbor_graph = "knn",
                                        cores = 2)
    head(cds_graph_test_results)
    deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE), ], q_value < 0.05))
    
    p1 <- plot_cells(cds,
                    genes = head(deg_ids),
                    show_trajectory_graph = FALSE,
                    label_cell_groups = FALSE,
                    label_leaves = FALSE)
    
    png(paste0(projDirS, 'Plots/Identify_genes_that_change_as_a_function_of_pseudotime.png'), width = 4000, height = 2500, res = 300)
    print(wrap_plots(p1))
    dev.off()
}
