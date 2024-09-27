library("stringr")
library('fgsea')
library('tidyr')

# # Set variable
# cluster_list = list('0.2','0.4', '0.8')
# known_cell_type_marker_file = '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_test1/MCP_PDGFB/mouse_markers.csv'#'/ahg/regevdata/projects/ICA_Lung/Nishant/LabCode/marker_file/mouse_celltypes_only_lymphoid_marker.csv'#'/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_test1/MCP_PDGFB/mouse_markers.csv'
# scrna_pipeline_dir = '/ahg/regevdata/projects/ICA_Lung/Nishant/LabCode/scrna_pipeline/'
# max_marker_for_Low_filter_file_known_marker = 2
# max_marker_for_high_filter_file_known_marker = 1
# max_marker_for_low_filter_file_no_cell_type_marker = 2
# max_marker_for_High_filter_file_no_cell_type_marker = 1
# Unknown_marker = TRUE ##  To include unknown marker
# force = FALSE ## To re run the complete pipeline
# outDir = denova_marker_path
# res_names = "orig.ident_harmony_snn_res."
# dotplot_1_group_by = 'sub_celltype_Compartment' ## It will create dotplot with markers, Select celltype column in seurat obj
# dotplot_2_group_by = 'orig.ident_harmony_snn_res.' ## It will create dotplot with markers, Select res column in seurat obj (only first common word)
# ## If seurat object is with other name
# # srt <- seurat object

# ## For GSEA analysis
# top_pathways = 10 ## Max Number of pathway select for dotplot 
# pvalAdjTrheshold = 0.05
# org='mouse' ## human or mouse
# ## If seurat object is with other name
# # srt <- seurat object
# source (paste0(scrna_pipeline_dir, 'Denova_study.R'))



source (paste0(usefull_function_path,'/useful_functions.R'), local = TRUE)

print('Start calculating marker for each cluter at diffrent res')
dir.create(paste0(outDir,'Plots/'))
for (cluster_res in cluster_list) {
    print(cluster_res)
    Idents(srt) <- paste0(res_names, cluster_res)

    if (file.exists(paste0(outDir, "srt",cluster_res,".rds")) & force){
        srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        saveRDS(srt.markers, file = paste0(outDir, "srt",cluster_res,".rds"))
        srtptwo <-readRDS(file = paste0(outDir, "srt",cluster_res,".rds"))
        Idents(srt) <- paste0(res_names, cluster_res)
        srtptwo%>%
        group_by(cluster) %>%
        top_n(n = 30, wt = avg_log2FC) -> top30
        write.csv(top30, paste0(outDir, "markers_file_cluster_", cluster_res, "_top30_gene.csv"), row.names=FALSE)
        pdf (paste0(outDir, 'Plots/newHEAT_MAP_RES_',cluster_res,'_legend_with_de_novo_marker.pdf'), width=90, height=100)
        print(DoHeatmap(srt, features = top30$gene)+ NoLegend()+ theme(text = element_text(size = 20)))
        dev.off()
        df_with_special_characters <- readRDS(paste0(outDir, "srt",cluster_res,".rds"))
        write.csv(df_with_special_characters, paste0(outDir, "markers_file_cluster_", cluster_res, ".csv"), row.names=FALSE)
    }
    else if (file.exists(paste0(outDir, "srt",cluster_res,".rds")) & !force) {
        print('file exit')
    } else{
        srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        saveRDS(srt.markers, file = paste0(outDir, "srt",cluster_res,".rds"))
        srtptwo <-readRDS(file = paste0(outDir, "srt",cluster_res,".rds"))
        Idents(srt) <- paste0(res_names, cluster_res)
        srtptwo%>%
        group_by(cluster) %>%
        top_n(n = 30, wt = avg_log2FC) -> top30
        write.csv(top30, paste0(outDir, "markers_file_cluster_", cluster_res, "_top30_gene.csv"), row.names=FALSE)
        pdf (paste0(outDir, 'Plots/newHEAT_MAP_RES_',cluster_res,'_legend_with_de_novo_marker.pdf'), width=90, height=100)
        print(DoHeatmap(srt, features = top30$gene)+ NoLegend()+ theme(text = element_text(size = 20)))
        dev.off()
        df_with_special_characters <- readRDS(paste0(outDir, "srt",cluster_res,".rds"))
        write.csv(df_with_special_characters, paste0(outDir, "markers_file_cluster_", cluster_res, ".csv"), row.names=FALSE)
    }

    # srtptwo <-readRDS(file = paste0(outDir, "srt",cluster_res,".rds"))
    # Idents(srt) <- paste0(res_names, cluster_res)
    # srtptwo%>%
    # group_by(cluster) %>%
    # top_n(n = 30, wt = avg_log2FC) -> top30
    # write.csv(top30, paste0(outDir, "markers_file_cluster_", cluster_res, "_top30_gene.csv"), row.names=FALSE)
    # pdf (paste0(outDir, 'Plots/newHEAT_MAP_RES_',cluster_res,'_legend_with_de_novo_marker.pdf'), width=90, height=100)
    # print(DoHeatmap(srt, features = top30$gene)+ NoLegend()+ theme(text = element_text(size = 20)))
    # dev.off()
    # df_with_special_characters <- readRDS(paste0(outDir, "srt",cluster_res,".rds"))
    # write.csv(df_with_special_characters, paste0(outDir, "markers_file_cluster_", cluster_res, ".csv"), row.names=FALSE)
}


print('Filter best marker for each cluster')
for (cluster_res in cluster_list) {
    Low_filter_for_known_marker = max_marker_for_Low_filter_file_known_marker
    high_filter_for_known_marker = max_marker_for_high_filter_file_known_marker
    low_filter_for_no_cell_type_marker = max_marker_for_low_filter_file_no_cell_type_marker
    High_filter_for_no_cell_type_marker = max_marker_for_High_filter_file_no_cell_type_marker

    res = paste0('Marker_select_at_res_',cluster_res)
    input_file = paste0('markers_file_cluster_',cluster_res,'.csv')
    input_file_dic = outDir

    system(paste0('python ',scrna_pipeline_dir, 'Immune_marker_select_pipeline.py ', input_file,' ',
    known_cell_type_marker_file,' ', res,' ',Low_filter_for_known_marker, ' ',high_filter_for_known_marker, ' ',
    low_filter_for_no_cell_type_marker, ' ',High_filter_for_no_cell_type_marker, ' ',input_file_dic), wait=TRUE)

}

Sys.sleep(5)

print('Create feature plot')
list.of.files <- list.files(outDir, "filter")
for (file_name in list.of.files) {
    print(file_name) 

    split_file_name = strsplit(file_name, '.csv', fixed=T)[[1]]
    print(split_file_name)

    markers.selected_df = read.csv(paste0(outDir, file_name))
    markers = markers.selected_df[!duplicated(markers.selected_df$gene),]
    markers = markers[markers$gene != '',]
    ### In denovo gene list 'DNMT1' is no available
    markers$gene[markers$gene == 'DNMT1'] = 'RNA_DNMT1'
    if (Unknown_marker == FALSE){
        markers = markers[markers$celltype != 'no_cell_type',]
    }

    message('Generate feature plots for the markers')
    if (length(markers$gene) > 300) {
        Font_size = 6
    } else if (length(markers$gene) > 200) {
        Font_size = 8
    } else if (length(markers$gene) > 100) {
        Font_size = 10        
    } else {
        Font_size = 12
    }
    # =======
    p = FeaturePlot(srt, features = markers$gene,
    ncol=4, order=F,col=viridis::turbo(100), combine = FALSE, pt.size = .001, reduction = reductionName)
    for(i in 1:length(p)) {
    p[[i]] = p[[i]] + theme_void() + NoAxes() + NoLegend() +
    ggtitle (paste (colnames(p[[i]]$data)[4],'\n', paste(markers[markers$gene == colnames(p[[i]]$data)[4],'celltype'], collapse = '_'))) +
    theme(plot.title = element_text(size = Font_size))
    }
    png(paste0(outDir,'Plots/',split_file_name,'_umaps.png'), 4000,4000, res=300)#2000,2000,Res=300(10000,10000, res=300)
    # >>>>>>> 4255fea734e96fe36a4ea0ba7ec587f5df0f0c5c
    print(wrap_plots(p))
    dev.off()

    ### Dotplot
    p1 <- DotPlot(object = srt, features = markers$gene, scale = T, group.by = dotplot_1_group_by) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
    scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

    res_no_split <- str_split(split_file_name, '_')[[1]][5]
    p2 <- DotPlot(object = srt, features = markers$gene, scale = T, group.by = paste0(dotplot_2_group_by,res_no_split)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
    scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

    if (length(markers$gene) <= 40) {
        dotplot_width = 15
    } else if (100 >= length(markers$gene)){
        dotplot_width = 35
    } else if(200 >= length(markers$gene)){
        dotplot_width = 45
    } else if(length(markers$gene) >= 200){
        dotplot_width = 65
    } else if(length(markers$gene) >= 300){
        dotplot_width = 85
    }

    if (length(table(srt@meta.data[,paste0(dotplot_2_group_by,res_no_split)])) <= 12){
        dotplot_hight = 8
    } else if (length(table(srt@meta.data[,paste0(dotplot_2_group_by,res_no_split)])) > 12){
        dotplot_hight = (length(table(srt@meta.data[,paste0(dotplot_2_group_by,res_no_split)])))/1.5
    }

    pdf (paste0(outDir,'Plots/',split_file_name,'_dotplot.pdf'), useDingbats = F, width = dotplot_width, height = dotplot_hight)
    print(p1 / p2)
    dev.off()

}

## Generate UMAP with addmodulescore of top 30 genes
for (cluster_no in cluster_list){
    top_30_genes <- read.csv(paste0(outDir, 'markers_file_cluster_',cluster_no,'_top30_gene.csv'))
    umap_p1_list <- list()
    for (i in 0:max(top_30_genes$cluster)){
        subset_gene <- subset(top_30_genes, cluster==i)
        modulesL = list (wgcna_modules=subset_gene$gene)
        message ('Run AddModuleScore')
        srt = AddModuleScore (srt, modulesL, name = names(modulesL))
        colnames (srt@meta.data)[colnames(srt@meta.data) %in% paste0(names(modulesL),seq(length(modulesL)))] = names (modulesL)
        meta_modules_names = names(modulesL)[names(modulesL) != 'grey'] # remove grey module (unassigned)
        umap_df = data.frame (srt[[reductionName]]@cell.embeddings, wgcna_modules = srt@meta.data[,meta_modules_names])


        umap_p1 <- lapply (meta_modules_names, function(x) ggplot(data = umap_df) + 
        geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
        scale_color_viridis (option='turbo')  +
        theme_bw() + ggtitle (paste0('cluster_',i)))

        umap_p1_list <- c(umap_p1_list, umap_p1)

        srt@meta.data[,'wgcna_modules'] <- NULL
    }
    if (max(top_30_genes$cluster) < 10){
        width_size = 8000
    } else if(max(top_30_genes$cluster) < 20) {
        width_size = 12000
    } else if(max(top_30_genes$cluster) < 30) {
        width_size = 16000
    } else if(max(top_30_genes$cluster) < 40) {
        width_size = 20000
    } else if(max(top_30_genes$cluster) < 50) {
        width_size = 24000
    } else if(max(top_30_genes$cluster) < 60) {
        width_size = 28000
    } else if(max(top_30_genes$cluster) >= 60) {
        width_size = 30000
    }
    png (paste0 (outDir,'Plots/addmodule_score_UMAP_with_top30_genes_at_res_',cluster_no,'.png'), width = width_size, height = 2000, pointsize=10, res = 300, type="cairo")
    print (wrap_plots (umap_p1_list, ncol = ifelse (length(umap_p1_list) > 8,ceiling(length(umap_p1_list)/2),length(umap_p1_list))))
    dev.off()
}


## GSEA analysis
file_output_path <- 'GSEA/'
plot_output_path <- paste0(file_output_path,'Plots/')
dir.create(paste0(outDir,file_output_path))
dir.create(paste0(outDir,plot_output_path))
do.fgsea = TRUE

for (cluster_no in cluster_list){
    print('Create GSEA dotplot')
    list.of.files <- list.files(outDir, paste0("srt",cluster_no,".rds"), full.names = FALSE)
    print(list.of.files)    
    for (srt_filenames in list.of.files){
        srt_find_marker <- readRDS(paste0(outDir,srt_filenames))
        if (do.fgsea)
        {
            if (org == 'human'){
                gmt_annotations = list(#'current_MSigDB_gene_sets_as_Gene_Symbols' = 'msigdb.v7.5.1.symbols.gmt', 
                                 #'chemical_and_genetic_perturbations_as_Gene_Symbols' = 'c2.cgp.v7.5.1.symbols.gmt',
                                #'all_canonical_pathways_as_Gene_Symbols' = 'c2.cp.v7.5.1.symbols.gmt',
                                'Reactome_gene_sets_as_Gene_Symbols' = 'c2.cp.reactome.v7.5.1.symbols.gmt',
                                #'all_regulatory_target_gene_sets_as_Gene_Symbols' = 'c3.all.v7.5.1.symbols.gmt',
                                 #'cancer_gene_neighborhoods_as_Gene_Symbols' = 'c4.cgn.v7.5.1.symbols.gmt',
                                 #'GO_biological_processes_as_Gene_Symbols' = 'c5.go.bp.v7.5.1.symbols.gmt',
                                 #'all_oncogenic_signature_gene_sets_as_Gene_Symbols' = 'c6.all.v7.5.1.symbols.gmt',
                                 #'all_immunologic_signature_gene_sets_as_Gene_Symbols' = 'c7.all.v7.5.1.symbols.gmt',
                                 #'all_cell_type_signature_gene_sets_as_Gene_Symbols' = 'c8.all.v7.5.1.symbols.gmt',
                                  'Go_gene' = 'c5.bp.v7.1.symbol.gmt',
                                  'hallmark_gene_sets_as_Gene_Symbols' = 'h.all.v7.5.1.symbols.gmt')
            } else if(org == 'mouse'){
                gmt_annotations = list('hallmark_gene_sets_as_Gene_Symbols' = 'mh.all.v2022.1.Mm.symbols.gmt',
                                       'Reactome_gene_sets_as_Gene_Symbols' = 'm2.cp.reactome.v2022.1.Mm.symbols.gmt',
                                       #'all_cell_type_signature_gene_sets_as_Gene_Symbols' = 'm8.all.v2022.1.Mm.symbols.gmt',
                                       #'Canonical_Pathways_gene_sets_derived_from_the_WikiPathways_database' = 'm2.cp.wikipathways.v2022.1.Mm.symbols.gmt',
                                       #'Canonical_Pathways_gene_sets_derived_from_the_BioCarta_pathway_database'='m2.all.v2022.1.Mm.symbols.gmt',
                                       'GO_Gene_Ontology_gene_sets' = 'm5.go.v2023.1.Mm.symbols.gmt'
                                       )
                }
            
            gmt_no = 1
            for (ann in gmt_annotations){
                    if (!file.exists (paste0(outDir, file_output_path,'fgsea_',ann,'_',cluster_no,'.rds')) | force){
                        fgseaResAll = list()
                        gmt.file = paste0 ('/ahg/regevdata/projects/ICA_Lung/Nishant/GSEA_DB/',org,'/',ann)
                        pathways = gmtPathways (gmt.file)
                        message (paste('Compute enrichment per cluster using annotation:', ann))
                        fgseaResCluster = list()
                        cluster_number <- 0
                        for (i in unique (srt_find_marker$cluster))
                            {
                            message (paste ('fGSEA running cluster',i))
                            tblCluster = srt_find_marker[srt_find_marker$cluster == i,]
                            tblCluster = na.omit (tblCluster)
                            fgsea_ranks = -log10 (tblCluster$p_val + 1e-300) * sign (tblCluster$avg_log2FC)
                            fgsea_ranks = setNames (fgsea_ranks, tblCluster$gene)
                            fgsea_ranks = fgsea_ranks[fgsea_ranks != 0]
                            fgseaRes = fgseaMultilevel (pathways, 
                                    fgsea_ranks, 
                                    minSize=15, 
                                    maxSize=1500, # changed this from 1500 to 1000 cause it generated an error
                                    BPPARAM = NULL)
                            fgseaResCol = collapsePathways (fgseaRes, stats = fgsea_ranks, pathway = pathways)
                            mainPathways = fgseaRes[fgseaRes$pathway %in% fgseaResCol$mainPathways]
                            fgseaResCluster[[i]] = mainPathways
                            mainPathways$cluster <- i
                            leadingEdge_gene <- mainPathways$leadingEdge
                            mainPathways$leadingEdge <- NULL
                            mainPathways$leadingEdge <- leadingEdge_gene
                            mainPathways1 <- apply(mainPathways,2,as.character)
                            if (length(dim(mainPathways1)) != 0){
                                if (cluster_number == 0){
                                    mainPathways1_data = mainPathways1
                                } else{
                                    mainPathways1_data <- rbind(mainPathways1_data,mainPathways1)
                                }
                            cluster_number = cluster_number+1
                            }
                            
                            }
                        write.csv(mainPathways1_data, paste0(outDir, file_output_path,'top_gene_with_pathway_',ann,'_',cluster_no,'.csv'), quote = F)
                        fgseaResAll[[ann]] = fgseaResCluster
                    saveRDS (fgseaResAll, paste0(outDir, file_output_path,'fgsea_',ann,'_',cluster_no,'.rds'))
                    } else {
                    message ('fGSEA object found!')
                    fgseaResAll = readRDS (paste0(outDir, file_output_path,'fgsea_',ann,'_',cluster_no,'.rds'))
                    }
                    
                    # Plot dotplot of fGSEA annotations per cluster 
                    fgseaResAll_dp = lapply (fgseaResAll, function(y) dotGSEA (y, padj_threshold = pvalAdjTrheshold, 
                        type = 'fgsea',top_pathways = top_pathways,
                        cluster_rows=T,
                        cluster_cols=T, 
                        title_name = paste0(names(gmt_annotations[gmt_no]), ' | Top pathway selected = ',top_pathways))
                    )
                    
                    if (length(unique(fgseaResAll_dp[[1]]$data$pathway)) < 20) {
                        fig_length <- 6
                    } else if (length(unique(fgseaResAll_dp[[1]]$data$pathway)) < 40) {
                        fig_length <- 10
                    } else if(length(unique(fgseaResAll_dp[[1]]$data$pathway)) < 70) {
                        fig_length <- 15
                    } else if(length(unique(fgseaResAll_dp[[1]]$data$pathway)) < 100) {
                        fig_length <- 20
                    } else if(length(unique(fgseaResAll_dp[[1]]$data$pathway)) < 140) {
                        fig_length <- 30
                    } else if (length(unique(fgseaResAll_dp[[1]]$data$pathway)) >= 140) {
                        fig_length <- 40
                    }
                        
                    pdf (paste0(outDir,plot_output_path,'fGSEA_',names(gmt_annotations[gmt_no]),'_',cluster_no,'_annotation_',names (fgseaResAll_dp)[1],'_dotplots.pdf'),15,fig_length)
                    print(fgseaResAll_dp[[1]])
                    dev.off()
                    gmt_no = gmt_no+1

            }
        }
    }

}

file_output_path <- 'GSEA_Enrichr/'
plot_output_path <- paste0(file_output_path,'Plots/')
dir.create(paste0(outDir,file_output_path))
dir.create(paste0(outDir,plot_output_path))
enricher_universe = rownames (srt)

for (cluster_no in cluster_list){
    print('Create GSEA dotplot')
    list.of.files <- list.files(outDir, paste0("srt",cluster_no,".rds"), full.names = FALSE)
    print(list.of.files)
    
    for (srt_filenames in list.of.files){
        srt_find_marker <- readRDS(paste0(outDir,srt_filenames))
            if (org == 'human'){
                gmt_annotations = list('hallmark_gene_sets_as_Gene_Symbols' = 'h.all.v7.5.1.symbols.gmt',
                                'Reactome_gene_sets_as_Gene_Symbols' = 'c2.cp.reactome.v7.5.1.symbols.gmt',
                                  'Go_gene' = 'c5.bp.v7.1.symbol.gmt'
                                  )
            } else if(org == 'mouse'){
                gmt_annotations = list('hallmark_gene_sets_as_Gene_Symbols' = 'mh.all.v2022.1.Mm.symbols.gmt',
                                       'Reactome_gene_sets_as_Gene_Symbols' = 'm2.cp.reactome.v2022.1.Mm.symbols.gmt',
                                       'GO_Gene_Ontology_gene_sets' = 'm5.go.v2023.1.Mm.symbols.gmt'
                                       )
                }
            # GSEA analysis on DEG per cluster
            EnrichRResAll = list()
            gmt_no = 1
            if (!file.exists (paste0(file_output_path, 'EnrichR_WGCNA_module_genes_',cluster_no,'.rds')) | force){
            for (ann in gmt_annotations){
                    #if (!file.exists (paste0(outDir, file_output_path,'fgsea_',ann,'_',cluster_no,'.rds')) | force){
                        #fgseaResAll = list()
                        gmt.file = paste0 ('/ahg/regevdata/projects/ICA_Lung/Nishant/GSEA_DB/',org,'/',ann)
                        pathways = read.gmt (gmt.file)
                        message (paste('Compute enrichment per cluster using annotation:', ann))
                        #fgseaResCluster = list()
                        EnrichRResCluster = list()
                        for (i in unique (srt_find_marker$cluster))
                            {
                            message (paste ('EnrichR running module',i))
                            #print(head(srt_find_marker))
                            sig_genes = srt_find_marker$gene[srt_find_marker$cluster == i]
                            egmt <- enricher (sig_genes, TERM2GENE=pathways, universe = enricher_universe)
                              if (length(egmt) == 0){
                                    message (paste ('Zoro genes match with Enricher')) 
                                } else{
                                  EnrichRResCluster[[i]] = egmt@result
                              }

                            }
                        EnrichRResAll[[ann]] = EnrichRResCluster
                        #}
                    message (paste ('Saving EnrichR ', ann))
                    saveRDS (EnrichRResAll, paste0(outDir, file_output_path, 'EnrichR_WGCNA_module_genes_',cluster_no,'.rds'))
                } 
                } else {
                message ('EnrichR object found!')
                EnrichRResAll = readRDS (paste0(outDir, file_output_path, 'EnrichR_WGCNA_module_genes_',cluster_no,'.rds'))
            }
                        
                    #pvalAdjTrheshold = 0.05
                    #top_pathways = 10
                    EnrichRes_dp = lapply (EnrichRResAll, function(x) dotGSEA (enrichmentsTest_list = x, type = 'enrich', 
                    padj_threshold = pvalAdjTrheshold, top_pathways= top_pathways))


                    for (i in seq_along(gmt_annotations))
                      {
                        if (i == 1){
                          width = 7
                          hight = 7
                        } else if (i==2){
                            width = 15
                            hight = 15
                        } else if (i==3){
                            width = 15
                            hight = 15
                        }

                      pdf (paste0(outDir, plot_output_path,'EnrichR_WGCNA_module_genes_',names(gmt_annotations[i]),'_',cluster_no,'_dotplots.pdf'),width, hight)
                      print (EnrichRes_dp[[i]])
                      dev.off()    

            }
    }
}




print('Complete!')