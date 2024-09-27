## Function available
### Calculate the cell density
print('Make "readme = TRUE" to see all function list')
readme <- TRUE
if (readme) {
        cat("To see the package readme make the 'helpme= TRUE' (eg. Cell_density_fun(helpme= TRUE)).\n")
        cat("This List of function available in this file.\n")
        cat("1. Default_fig() : Generate Default figures\n")
        cat("2. Enrichr_fun() : Pathway enrichment analysis with Enrichr\n")
        cat("3. feature_expression_analysis_function() : features expression analysis function - dotplot, feature plot, violin plot\n")
        cat("4. Cell_density_fun(): Calculate the cell density and make UMAP\n")
        cat("5. Valcano_plot(): It can calculate the DE gene between the given column ana make the valcano plot\n")
        return(invisible())
    }


############################################################################
############################################################################
### Check packages

pck_check <- function(required_packages){
    lapply(required_packages, require, character.only = TRUE)
    # Check which packages are not installed
    missing_packages <- setdiff(required_packages, rownames(installed.packages()))
    
    # Load missing packages if any
    if (length(missing_packages) > 0) {
        message("The following packages are missing:", paste(missing_packages, collapse = ", "))
        install_choice <- readline("Do you want to install these packages? (y/n): ")
        
        if (tolower(install_choice) == "y") {
        install.packages(missing_packages)
        } else {
        stop("Required packages are not installed. Function cannot proceed.")
        }
    }

    for (pkg in missing_packages) {
        library(pkg, character.only = TRUE)
    }
}



############################################################################
############################################################################
### usefull function
cellComp = function (
                    seurat_obj = NULL,
                    metaGroups = NULL, # vector of at least 3 metaGroups e.g. c('orig.ident','celltypes','celltypes'),
                    plot_as = c('box','bar'),
                    pal = NULL,
                    prop = TRUE,
                    ptable_factor = 1, # specify which column of the data.frame or seurat object metadata should be used to compute proportions
                    facet_ncol = 20,
                    facet_scales = 'free',
                    subset_prop = NULL, # subset prop table by any group in any column
                    pair_com = NULL,
                    Pvalue_cal = NULL,
                    Pvalue_method = 't_test',
                    hide_ns = NULL,
                    Pvalue_label =NULL,
                    p.adjust.method = NULL
)
    
{

if (is.data.frame (seurat_obj))
{
meta_groups_df = seurat_obj[,metaGroups]
} else {
meta_groups_df = seurat_obj@meta.data[,metaGroups]
}
if(is.null(pal)) pal = rainbow (length(unique(meta_groups_df[,2])))
if (prop)
{
ccomp_df = as.data.frame (prop.table(table (meta_groups_df),ptable_factor))
} else {
ccomp_df = as.data.frame (table (meta_groups_df))
}
ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
if (!is.null (subset_prop))
{
subset_col = unlist(sapply (seq(ncol(ccomp_df)), function(x) if(any(ccomp_df[,x] %in% subset_prop)) colnames(ccomp_df)[x]))
ccomp_df = ccomp_df[ccomp_df[,subset_col] %in% subset_prop,]
}
if (plot_as == 'box')
    {

    colnames (ccomp_df) = c(paste0('Var_',seq_along(metaGroups)), 'proportion')


    box_p = ggplot (ccomp_df, aes (x= Var_2, y= proportion)) +
    geom_boxplot(aes (fill= Var_3)) +
    theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual (values= pal) #+
    #geom_pwc(method = "t_test", label = "{p.format}{p.signif}", hide.ns = T, p.adjust.method = "none")

    if (Pvalue_cal){
        box_p <- box_p + geom_pwc(method = Pvalue_method, label = Pvalue_label, hide.ns = hide_ns, p.adjust.method = p.adjust.method)
    }


    if (length(metaGroups) > 3) box_p = box_p + facet_wrap (~Var_4, scales=facet_scales, ncol=facet_ncol)
}
if (plot_as == 'bar')
    {
    colnames (ccomp_df) = c(paste0('Var_',seq_along(metaGroups)), 'proportion')
    box_p = ggplot (ccomp_df, aes (x= Var_1, y= proportion)) +
    geom_bar(position="stack", stat="identity", aes(fill= Var_2)) +
    theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual (values= pal)
    if (length(metaGroups) == 3) box_p = box_p + facet_wrap (~Var_3, scales=facet_scales, ncol=facet_ncol)
    }
return (box_p)
}

############################################################################
############################################################################
### Default figure

Default_fig <- function(srt =NULL, fig_no= 'ALL', out_path = NULL, sampleid_col = NULL, 
                        Genotype_col = NULL, percent.mt = NULL, res = NULL, org = NULL, 
                        ccs_ccg2m_threshold=NULL, col_list_for_UMAP = NULL, 
                        usefull_function_path = NULL, cluster_list =list('0.2', '0.4','0.6', '0.8', '1.2', '2', '3'), 
                        known_cell_type_marker_file = NULL,
                        scrna_pipeline_dir = NULL, force = FALSE, res_names_col =NULL, 
                        dotplot_1_group_by_col = NULL, pvalAdjTrheshold = 0.05,
                        top_pathways = 10, helpme = FALSE){

        

    # Check if helpme parameter is TRUE
    if (helpme) {
        cat("This function generate the following Default fig.\n")
        cat("\n\nHere is the list figure generated by this function.\n")
        cat("1. nfeat ncount and mpercentage vln plot at genotype level\n")
        cat("2. nfeat ncount and mpercentage vln plot at sample level\n")
        cat("3. nCounts and nfeature umaps\n")
        cat("4. Percentage of mitochondria umap\n")
        cat("5. De-novo cluster umap at diffrent res\n")
        cat("6. Cell cycle analysis\n")
        cat("7. Generate UMAP at diffrent metadata\n")
        cat("8. Denovo maker identification - This function required Denovo_study_ns.R and  Immune_marker_select_pipeline.py script in scrna_pipeline folder\n")

        cat("\n\n\n")

        cat("Required input parameters:\n")
        cat("For all figure:\n")
        cat("- srt: seurat object.\n")
        cat("- fig_no: c(1, 2, 5) to generate specific fig or 'ALL' for generating the all figure\n")
        cat("- out_path: Output path for figutre (eg. paste(projDir/Plots/))\n")
        cat("Parameter required For fun-1:\n")
        cat("- Genotype_col: Genotype column name\n")
        cat("Parameter required For fun-2:\n")
        cat("- sampleid_col: sampleid column name\n")
        cat("Parameter required For fun-4:\n")
        cat("- percent.mt: percent.mt column name\n")      
        cat("Parameter required For fun-5:\n")  
        cat("- res: Resolution list (eg. c(0.2, 0.4))\n")
        cat("Parameter required For fun-6:\n")
        cat("- org: 'mouse' or 'human'\n")
        cat("- ccs_ccg2m_threshold: c(0.10, 0.10)\n")
        cat("Parameter required For fun-7:\n")
        cat("- col_list_for_UMAP: c('celltype', 'sub_celltype')\n")
        cat("Parameter required For fun-8:\n")
        cat("- usefull_function_path: usefull_function script path\n")
        cat("- cluster_list: cluster res list\n")
        cat("- known_cell_type_marker_file: known cell_type marker file path (contain gene and celltype colnames) \n")
        cat("- scrna_pipeline_dir: scrna pipeline dir path to read the Denovo_study_ns.R and  Immune_marker_select_pipeline.py \n")
        cat("- force: TRUE and FALSE\n")
        cat("- res_names_col: res Col name in seurat obejct without res velue for read all res column (eg sampleID_harmony_snn_res.) \n")
        cat("- dotplot_1_group_by_col: dotplot 1 y axis (eg . 'celltype')\n")
        cat("- pvalAdjTrheshold: pvalAdjTrheshold value for pathway enrichment (eg . 0.05)\n")
        cat("- top_pathways: top_pathways value (eg . 10)\n")

        cat("eg: Default_fig(srt =srt, fig_no= 8, out_path = out_path, sampleid_col = 'sampleID', 
                        Genotype_col = 'groupID', percent.mt = 'percent.mt', res = res, org = 'mouse', 
                        ccs_ccg2m_threshold=c(0.10, 0.10), col_list_for_UMAP = c('celltype'), 
                        usefull_function_path = usefull_function_path, cluster_list =list('0.2', '0.4','0.6', '0.8', '1.2', '2'), 
                        known_cell_type_marker_file = known_cell_type_marker_file,
                        scrna_pipeline_dir = scrna_pipeline_dir, force = FALSE, res_names_col ='sampleID_harmony_snn_res.', 
                        dotplot_1_group_by_col = c('celltype'), pvalAdjTrheshold = 0.05,
                        top_pathways = 10)\n")

        
        return(invisible())
    }

    # List of required packages
    packages = c('ggplot2','Seurat','dplyr','Matrix','ggpubr',
        'gplots','patchwork','ComplexHeatmap',
        'RColorBrewer','ggrepel','fgsea',
        'tidyr','gdata','GSA',
        'harmony','scater',
        'destiny','plotly',
        'GO.db', 'org.Mm.eg.db','org.Hs.eg.db',
        'viridis',
        'BiocParallel',
        'ggridges',
        'clusterProfiler','igraph','tidyverse', 'paletteer', 'rstatix', 'patchwork'
        )
    pck_check(packages)


    if (any(fig_no %in% c(1)) | any(fig_no == 'ALL')) {
    ## to generate nfeat ncount and mpercentage vln plot
    message ('nfeat ncount and mpercentage vln plot')

    Idents (srt ) = srt@meta.data[,Genotype_col]
    ccomp_df = as.data.frame (table (srt@meta.data[,Genotype_col]))
    cc_p2 = ggplot (ccomp_df, aes (x= Var1, y= Freq)) +
        geom_bar (position="stack", stat="identity") +
        theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=0.5, size = 10)) + 
        theme (axis.text = element_text (size = 12)) +
        ggtitle (paste('Total cells',ncol(srt)))+
        theme(plot.title = element_text(hjust = 0.5, size = 15, face="bold"))

    vln_p1 = VlnPlot (srt, features = "nFeature_RNA", , group.by = Genotype_col,pt.size = 0, ncol = 1)+
        theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=0.5, size = 10)) +
        theme(legend.position="none")
    vln_p2 = VlnPlot (srt, features = "nCount_RNA", group.by = Genotype_col,pt.size = 0, ncol = 1)+
        theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=0.5, size = 10)) +
        theme(legend.position="none")
    vln_p3 = VlnPlot (srt, features = "percent.mt", group.by = Genotype_col,pt.size = 0, ncol = 1)+
        theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=0.5, size = 10)) +
        theme(legend.position="none")

    png (paste0(out_path, "QC_nFeat_nCount_m.percent_vlnPlot_at_Genotype_level.png"), 4000, 1500, res=300)
    print ((cc_p2 + vln_p1 + vln_p2 + vln_p3) + plot_layout(widths=c(1,1,1,1)))
    dev.off()
    }

    if (any(fig_no %in% c(2)) | any(fig_no == 'ALL')) {
    ## to generate nfeat ncount and mpercentage vln plot
    message ('nfeat ncount and mpercentage vln plot')
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

    png (paste0(out_path, "QC_nFeat_nCount_m.percent_vlnPlot_at_sample_level.png"), 5000, 1500, res=300)
    print ((cc_p2 + vln_p1 + vln_p2 + vln_p3) + plot_layout(widths=c(1,1,1,1)))
    dev.off()
    }



    if (any(fig_no %in% c(3)) | any(fig_no == 'ALL')) {
        ## ncounts and nfeature umaps
        message ('nCounts and nfeature umaps')
        metaGroupNames = c(Genotype_col,sampleid_col)
        ####
        humap_p = lapply (metaGroupNames, function(y) DimPlot (object = srt, reduction = reductionName, pt.size = 1, group.by = y))

        # Color UMAP by number of detected genes and color from 5% to 95% of counts
        umap_df = data.frame (srt[[reductionName]]@cell.embeddings, nCount_RNA = log10(srt$nCount_RNA+1), nFeature_RNA = log10(srt$nFeature_RNA+1))
        umap_p1 = ggplot(data = umap_df) + 
        geom_point(mapping = aes_string(x = colnames(umap_df)[1], y=colnames(umap_df)[2], color = colnames(umap_df)[3]), size = .01) + 
        scale_colour_gradientn (colours = rainbow(7)) +
        theme_bw() +
        theme(plot.background = element_blank())+
        theme(legend.position="bottom")+
        theme(axis.text.x=element_text(size=15))+
        theme(axis.text.y=element_text(size=15))

        umap_p2 = ggplot(data = umap_df) + 
        geom_point(mapping = aes_string(x = colnames(umap_df)[1], y=colnames(umap_df)[2], color = colnames(umap_df)[4]), size = .01) + 
        scale_colour_gradientn (colours = rainbow(7)) +
        theme_bw() +
        theme(plot.background = element_blank())+
        theme(legend.position="bottom")+
        theme(axis.text.x=element_text(size=15))+
        theme(axis.text.y=element_text(size=15))
        png (paste0(out_path, 'Merged_nCounts_nFeatures_umap.png'), width = 5000, height = 3000, pointsize=10, res = 300, type="cairo")
        print (wrap_plots (c(humap_p, list(umap_p1, umap_p2))))
        dev.off()
    }


    if (any(fig_no %in% c(4)) | any(fig_no == 'ALL')) {
        ## percentage mitrochondria umap
        message ('Percentage mitochondria umap')
        # Color UMAP by number of detected genes and color from 5% to 95% of counts
        umap_df = data.frame (srt[[reductionName]]@cell.embeddings, percent.mt = srt$percent.mt)
        umap_p1 = ggplot(data = umap_df) + 
        geom_point(mapping = aes_string(x = colnames(umap_df)[1], y=colnames(umap_df)[2], color = colnames(umap_df)[3]), size = .01) + 
        scale_colour_gradientn (colours = rainbow(7)) +
        theme_bw() +
        theme(
        plot.background = element_blank()
        )

        png (paste0(out_path, 'percent.mt.png'), width = 2000, height = 1500, pointsize=10, res = 300, type="cairo")
        print (wrap_plots (umap_p1))
        dev.off()
    }



    if (any(fig_no %in% c(5)) | any(fig_no == 'ALL')) {
        ### Denovo cluster umap at diffrent res
        message ('Denovo cluster umap at diffrent res')
        clust_p1 = list()
        for (i in seq_along(res))
            {
            clust_p1[[i]] = DimPlot (srt, pt.size = 1, label = T, group.by= paste0(reductionGraph,'_res.',res[i]), reduction = reductionName) + NoLegend()
            }
        png (paste0(out_path, 'denovo_clusters_',reductionSave,'_umaps.png'), 3000, 4000, pointsize=10, res = 300, type="cairo")
        print (wrap_plots (clust_p1, ncol=2))
        dev.off()
    }


    if (any(fig_no %in% c(6)) | any(fig_no == 'ALL')) {
        ### Cell cycle analysis
        message ('Cell cycle analysis')
        if (org == 'mouse') {
            m.g2m.genes = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/cellcycle_mouse.g2m.genes.rds')
            m.s.genes = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/cellcycle_mouse.s.genes.rds')
            s.genes <- m.s.genes
            g2m.genes <- m.g2m.genes
        } else if (org == 'human'){
            s.genes <- cc.genes$s.genes
            g2m.genes <- cc.genes$g2m.genes
        }
        
        ## Calculate CellCycleScoring score
        srt <- CellCycleScoring(srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
        
        metaGroupNames = c('Phase')
        umap <- DimPlot (object = srt, reduction = reductionName, pt.size = 1, cols =paletteer_c("grDevices::rainbow", dim(table(srt@meta.data[,metaGroupNames]))), group.by = metaGroupNames, label = F) +theme(legend.position="bottom")
        png (paste0(out_path, 'Cell_cycle_phase_umap.png'), width = 2200, height = 2200, pointsize=10, res = 300, type="cairo")
        print (wrap_plots (umap))
        dev.off()


        # Color UMAP by number of detected genes and color from 5% to 95% of counts
        umap_df = data.frame (srt[[reductionName]]@cell.embeddings, S.Score = srt$S.Score, G2M.Score = srt$G2M.Score)
        umap_p1 = ggplot(data = umap_df) + 
        geom_point(mapping = aes_string(x = colnames(umap_df)[1], y=colnames(umap_df)[2], color = colnames(umap_df)[3]), size = .01) + 
        scale_colour_gradientn (colours = rainbow(10)) +
        theme_bw() +
        theme(
        plot.background = element_blank()
        )
        umap_p2 = ggplot(data = umap_df) + 
        geom_point(mapping = aes_string(x = colnames(umap_df)[1], y=colnames(umap_df)[2], color = colnames(umap_df)[4]), size = .01) + 
        scale_colour_gradientn (colours = rainbow(10)) +
        theme_bw() +
        theme(
        plot.background = element_blank()
        )
        png (paste0(out_path,'Cellcycle_s_g2_m_score.png'), width = 2200, height = 1000, pointsize=10, res = 300, type="cairo")
        print (wrap_plots (umap_p1, umap_p2))
        dev.off()



        ### Thrushold value
        ccs_threshold = ccs_ccg2m_threshold[1]
        ccg2m_threshold = ccs_ccg2m_threshold[2]

        # histogram cell cycle score distribution
        ccc_sc = ggplot(as.data.frame(srt$S.Score), aes(x=srt$S.Score)) +
        geom_histogram() +
        geom_vline(xintercept=ccs_threshold, color = "red",size=2, linetype = "dotted") + theme_bw()
        png (paste0(out_path,'cell_cycle_histogram_s_score.png'), width = 2000, height = 1800, pointsize=10, res = 300, type="cairo")
        print (wrap_plots (ccc_sc))
        dev.off()

        ccc_gc = ggplot(as.data.frame(srt$G2M.Score), aes(x=srt$G2M.Score)) +
        geom_histogram() +
        geom_vline(xintercept=ccg2m_threshold, color = "red", size=2, linetype = "dotted") + theme_bw()
        png (paste0(out_path,'cell_cycle_histogram_g2m_score.png'), width = 2000, height = 1800, pointsize=10, res = 300, type="cairo")
        print (wrap_plots (ccc_gc))
        dev.off()

        ## Classify the cells into cycling and non-cycling
        srt$Phase2 = ifelse (srt$S.Score > ccs_threshold | srt$G2M.Score > ccg2m_threshold, 'cycling', 'non-cycle')

        ### Create UMAP
        ## Cycling color
        non_cycling_color <- "#CCCCCC"  # Light Gray for non-cycling cells
        cycling_color <- "#66c2a5"  # Green for cycling cells
        cc.color <- c(cycling_color, non_cycling_color)
        metaGroupNames = c('Phase2')
        umap <- DimPlot (object = srt, reduction = reductionName, pt.size = .1, cols =cc.color, group.by = metaGroupNames) #+theme(legend.position="bottom")
        png (paste0(out_path,'classify_cell_cycle_umap.png'), width = 2000, height = 1800, pointsize=10, res = 300, type="cairo")
        print (wrap_plots (umap))
        dev.off()

    }

    if (any(fig_no %in% c(7)) | any(fig_no == 'ALL')) {
        ### Denovo marker study
        message ('Make UMAP at diffrent metadata')
        for (colname in col_list_for_UMAP){
            metaGroupNames = c(colname)
            umap <- DimPlot (object = srt, reduction = reductionName, pt.size = 1, label = TRUE, cols =paletteer_c("grDevices::rainbow", dim(table(srt@meta.data[,metaGroupNames]))), group.by = metaGroupNames) #+theme(legend.position="bottom")
            png (paste0(projDir,'Plots/',colname,'_umap.png'), width = 3500, height = 2500, pointsize=10, res = 300, type="cairo")
            print (wrap_plots (umap))
            dev.off()
        }
    }


    if (any(fig_no %in% c(8)) | any(fig_no == 'ALL')) {
        ### Denovo marker study
        message ('Denovo markers study at different res')
        denova_marker_path <- paste0(out_path,'denova_marker_study/')
        dir.create(denova_marker_path)
        usefull_function_path <- usefull_function_path

        cluster_list = cluster_list## Cluster resolution list
        known_cell_type_marker_file = known_cell_type_marker_file# Know marker file, Column name should be with 'gene' and 'celltype'
        scrna_pipeline_dir = scrna_pipeline_dir
        max_marker_for_Low_filter_file_known_marker = 2
        max_marker_for_high_filter_file_known_marker = 1
        max_marker_for_low_filter_file_no_cell_type_marker = 2
        max_marker_for_High_filter_file_no_cell_type_marker = 1
        Unknown_marker = TRUE ##  To include unknown marker
        force = force ## To re run the complete pipeline
        outDir = denova_marker_path
        res_names = res_names_col
        dotplot_1_group_by = dotplot_1_group_by_col
        dotplot_2_group_by = res_names_col
        # reductionName = 
        ### If seurat object is with other name
        # srt <- seurat object
        ## Parameter for GSEA analysis
        pvalAdjTrheshold = 0.05
        top_pathways <- 10
        source (paste0(scrna_pipeline_dir, 'Denovo_study_ns.R'), local = TRUE)
    }

return(srt)
}



############################################################################
############################################################################
### Enrichr pathway analysis - Pathway enrichment analysis with Enrichr
Enrichr_fun <- function(srt = NULL, inputDir = NULL, org = NULL, 
                        meta_col_names = NULL, Findallmarker_file = NULL, force = TRUE,
                        pvalAdjTrheshold = 0.05, top_pathways = 10,
                        helpme = FALSE){

    # Check if helpme parameter is TRUE
    if (helpme) {
        cat("This function generate Enrichr pathway analysis dotplot\n")
        cat("Required input parameters:\n")
        cat("- srt: seurat object.\n")
        cat("- inputDir: Input dir path (eg. projDir).\n")
        cat("- org: 'mouse' or 'human'.\n")
        cat("- meta_col_names: meta column name select from seurat object for looking enrichr pathway analysis (eg. sub_celltype).\n")
        cat("- Findallmarker_file: provied Findallmarker_file names in list with complete path, if doest exit than provide 'NULL'.\n")
        cat("- pvalAdjTrheshold: pvalAdjTrheshold value for selecting pathway.- default is 0.05 \n")
        cat("- top_pathways: select top pathway for each suntype - default is 10\n")
        cat(" eg. Enrichr_fun(srt = srt, inputDir = inputDir, org = org,
            meta_col_names = meta_col_names, Findallmarker_file = Findallmarker_file, force = TRUE,
            pvalAdjTrheshold = 0.05, top_pathways = 10)")
        return(invisible())
    }
    
    file_output_path <- paste0(inputDir,'/GSEA_Enrichr/')
    plot_output_path <- paste0(file_output_path,'Plots/')
    dir.create(paste0(file_output_path))
    dir.create(paste0(plot_output_path))
    enricher_universe = rownames (srt)
    
    if (is.null(Findallmarker_file)){
        for (meta_col_name in meta_col_names){
            if (!file.exists (paste0(file_output_path, 'srt_top_markers_between_',meta_col_names,'.rds')) | force){
                Idents(srt) <- meta_col_names
                # Run differential expression analysis using FindAllMarkers
                srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

                # Save the marker file
                Findallmarker_file <- paste0(file_output_path, "srt_top_markers_between_",meta_col_names,".rds")
                saveRDS(srt.markers, file = Findallmarker_file)
            }
        }
        list.of.files <- list.files(file_output_path, pattern = paste0("srt_top_markers_between_", meta_col_names, "\\.rds$"), full.names = FALSE)
    } else {
        list.of.files <- Findallmarker_file
    }
    
    for (srt_filenames in list.of.files){
        srt_find_marker <- readRDS(paste0(file_output_path,srt_filenames))
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
            if (!file.exists (paste0(file_output_path, 'EnrichR_at_',meta_col_names,'.rds')) | force){
            for (ann in gmt_annotations){
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
                    saveRDS (EnrichRResAll, paste0(file_output_path, 'EnrichR_at_',meta_col_names,'.rds'))
                } 
                } else {
                message ('EnrichR object found!')
                EnrichRResAll = readRDS (paste0(file_output_path, 'EnrichR_at_',meta_col_names,'.rds'))
            }

                    #pvalAdjTrheshold = 0.05
                    #top_pathways = 10
                    EnrichRes_dp = lapply (EnrichRResAll, function(x) dotGSEA (enrichmentsTest_list = x, type = 'enrich', padj_threshold = pvalAdjTrheshold, top_pathways= top_pathways))

                    for (i in seq_along(gmt_annotations))
                      {
                        if (i == 1){
                          width = 7
                          hight = 7
                          width1 = 2000
                          hight1 = 2000
                        } else if (i==2){
                            width = 15
                            hight = 15
                            width1 = 3500
                            hight1 = 3500
                        } else if (i==3){
                            width = 15
                            hight = 15
                            width1 = 3000
                            hight1 = 3000
                        }

                      pdf (paste0(plot_output_path,'EnrichR_at_',names(gmt_annotations[i]),'_',meta_col_names,'_dotplots.pdf'),width, hight)
                      print (EnrichRes_dp[[i]])
                      dev.off()
                        
                      png (paste0(plot_output_path,'EnrichR_at_',names(gmt_annotations[i]),'_',meta_col_names,'_dotplots.png'),width1, hight1, res=300)
                      print (EnrichRes_dp[[i]])
                      dev.off()
                        
                    }
        }
    
}




############################################################################
############################################################################
### features expression analysis function - dotplot, feature plot, violin plot   

# gene = c('Cd274', 'Marco')
# gene_list <- list('test' = gene)
# srt <- srt
# genotype_colnames <- 'Genotype_DM' ## If want to see dot plot at sample level thant replace with sampleid
# celltype_colnames <- 'sub_celltype'
# reductionName <- 'sampleID_harmony_umap' ## Reduction name in seurat object
# pt.size = 0.1 ## Feature plot dot size
# out_path <- paste0(projDir, 'Plots/')

feature_expression_analysis_function <- function(srt = NULL, gene_list=NULL, genotype_colnames =NULL, 
                                                 celltype_colnames=NULL, reductionName = 'umap', pt.size = 0.1,
                                                out_path = NULL, helpme = FALSE) {

    # Check if helpme parameter is TRUE
    if (helpme) {
        cat("This function do the features expression analysis and generate the dotplot, feature plot, violin plot\n")
        cat("Required input parameters:\n")
        cat("- srt: seurat object.\n")
        cat("- gene_list: genes list (eg. list('test' = gene)) .\n")
        cat("- genotype_colnames: Genotype column name in seurat object .\n")
        cat("- celltype_colnames: Celltype column name in seurat object .\n")
        cat("- reductionName : reductionName name in seurat object .\n")
        cat("- pt.size: Feature plot dot size eg. 0.1\n")
        cat("- out_path: Output path for plots .\n")
        cat(" eg. feature_expression_analysis_function(srt = srt, gene_list = gene_list, genotype_colnames = genotype_colnames, 
                                    celltype_colnames = celltype_colnames, reductionName = reductionName,
                                    out_path = out_path)")
        return(invisible())
    } 

    gene_plot <- paste0(out_path, names(gene_list), '/')
    dir.create(gene_plot)
    
    ## Generate dotplot between celltype and genotype and feature plot
    for (gene_name in names(gene_list)){
        #print(gene_name)
        #print(gene_list[gene_name])
        message(gene_name)
        genes = gene_list[[gene_name]]
        names <- gene_name
        # Create lists for dot plots and feature plots
        plot_list = list()
        dotp = list()
        for (g in genes)
          {
            message(g)
          dotp[[g]] <- geneDot (
          mat_norm = srt@assays$RNA@data,
          mat_counts = srt@assays$RNA@counts,
          #gene = top_tfs2,
          gene = g,
          x = srt@meta.data[,genotype_colnames], # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
          y = srt@meta.data[,celltype_colnames],
          min_expression = 0,
          x_name = genotype_colnames,
          y_name = celltype_colnames,
          facet_ncol = 6,
          lim_expression = NULL,
          scale.data=FALSE,
          plotcol = rev(brewer.pal(11,"Spectral"))) #rev(viridis::mako(100)))


            # FeaturePlots of Il1a and Il1b
            featp = FeaturePlot (srt, features = g,
                               ncol=1, cols = viridis::turbo(100), combine=F, pt.size = 0.1, reduction = reductionName)
            for(i in 1:length(featp)) {
                featp[[i]] = featp[[i]] + NoLegend() + NoAxes()
                }
            plot_list[[g]] <- dotp[[g]] + featp
          }
        
        if (length(gene_list[[1]])  > 1){
            png (paste0(gene_plot, 'Dotplot_',names,'.png'), width=5000, height=1500*length(genes)/2, res=300)
            print(wrap_plots(plot_list, ncol=2, nrow=length(genes)/2))
            dev.off()
        } else {
            png (paste0(gene_plot, 'Dotplot_',names,'.png'), width=2500, height=1500, res=300)
            print(wrap_plots(plot_list, ncol=1, nrow=1))
            dev.off()
        }

    }
    
    ## Generate dotplot at celltype level
    genes <- unique(gene_list[[1]])
    if (length(genes) < 3){
        width = 7
    } else{
        width = 5+length(genes) 
    }
    if (length(unique(srt@meta.data[,celltype_colnames])) < 5){
        height = 4
    } else{
        height = 1+length(unique(srt@meta.data[,celltype_colnames]))
    }
    
    p1 <- DotPlot(object = srt, features = genes, scale = T, group.by = celltype_colnames) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
    scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

    pdf (paste0(gene_plot, names(gene_list),'_markers_dotplot_',celltype_colnames,'.pdf'), useDingbats = F, width = 5+length(genes), height = 2+length(celltype_colnames))
    print(wrap_plots(p1))
    dev.off()
    
       
    ## ## Generate dotplot at celltype level and splitted at genotype level
    Idents(srt) <- genotype_colnames
    dotplot_list = list()
    for (genoname in names(table(srt@meta.data[,genotype_colnames]))){        
        condition <- srt@meta.data[,genotype_colnames] == genoname
        srt2 <- srt[,condition]
        
        p1 <- DotPlot(object = srt2, features = genes, scale = T, group.by = celltype_colnames) +
        theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
        scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
        ggtitle(genoname)
        dotplot_list[[genoname]] <- p1
    }

    pdf (paste0(gene_plot, names(gene_list),'_markers_dotplot_at_',celltype_colnames,'_splited_at_',genotype_colnames,'.pdf'), useDingbats = F, width = 13, height = 8)
    print(wrap_plots(dotplot_list))
    dev.off()
    
    
    ## ## Generate heatmap at celltype level
    dh <- DoHeatmap(srt, features = genes, group.by=celltype_colnames, size = 3, angle = 30)
    png (paste0(gene_plot, 'Heatmap_genes_display_at_',celltype_colnames,'.png'), width = 4000, height = 4000, pointsize=10, res = 300, type="cairo")
    print (dh)
    dev.off()
    
    
    ## ## Generate Violin plot at celltype level
    Idents(srt) <- celltype_colnames
    vlpt <- VlnPlot(srt, features = gene)
    png (paste0(gene_plot, 'Violin_plot_display_at_',celltype_colnames,'.png'), width = 4000, height = 2000, pointsize=10, res = 300, type="cairo")
    print (wrap_plots(vlpt))
    dev.off()
    
}

## Genedot function required in feature_expression_analysis_function function
geneDot <- function (mat_norm = NULL, mat_counts = NULL, gene = NULL, x = NULL, 
    y = NULL, z = "none", x_name = "genes", y_name = "clusters", 
    min_expression = 0, facet_ncol = 5, lim_expression = NULL, 
    scale.data = TRUE, plotcol = rev(viridis::rocket(100)), ...) 
{
    gene = gene[gene %in% rownames(mat_norm)]
    if (length(gene) > 1) 
        x = "gene"
    xyz_facs = paste(x, y, z, sep = "|||")
    mat_counts[is.na(mat_counts)] = 0
    isExp_mat = mat_counts[gene, , drop = F] > min_expression
    exp_mat = mat_norm[gene, , drop = F]
    pch_expressed = unlist(as.data.frame(sapply(unique(xyz_facs), 
        function(x) Matrix::rowSums(isExp_mat[, xyz_facs == x, 
            drop = F])/ncol(isExp_mat[, xyz_facs == x, drop = F]))), 
        use.names = F)
    exp_mat[!isExp_mat > 0] = NA
    avg_expressed = unlist(as.data.frame(sapply(unique(xyz_facs), 
        function(x) Matrix::rowMeans(exp_mat[, xyz_facs == x, 
            drop = F], na.rm = T))), use.names = F)
    avg_expressed[is.na(avg_expressed)] = 0
    if (length(gene) > 1) {
        x_axis = rep(gene, length(unique(xyz_facs)))
    }
    else {
        x_axis = sapply(unique(xyz_facs), function(x) unlist(strsplit(x, 
            "\\|\\|\\|"))[1])
    }
    dot.df = data.frame(x_axis = x_axis, y_axis = rep(sapply(unique(xyz_facs), 
        function(x) unlist(strsplit(x, "\\|\\|\\|"))[2]), each = length(gene)), 
        z_axis = rep(sapply(unique(xyz_facs), function(x) unlist(strsplit(x, 
            "\\|\\|\\|"))[3]), each = length(gene)), expression = avg_expressed, 
        percent_expressed = pch_expressed)
    dot.df$percent_expressed = dot.df$percent_expressed * 100
    dot.df2 = dot.df
    if (scale.data == TRUE) 
        dot.df2 = transform(dot.df2, expression = ave(expression, 
            x_axis, FUN = function(x) scale(x, scale = T, center = T)))
    p = ggplot(data = dot.df2, aes(x = x_axis, y = y_axis)) + 
        geom_point(shape = 21, aes(fill = expression, size = percent_expressed), 
            colour = "black") + scale_fill_gradientn(colours = plotcol) + 
        labs(x = x_name, y = y_name, title = ifelse(length(gene) > 
            1, "", gene), subtitle = paste("Min expression >", 
            min_expression)) + scale_shape(solid = FALSE) + theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
            hjust = 1))
    if (!z[1] == "none") 
        p = p + facet_wrap(~z_axis, drop = TRUE, scales = "free_x", 
            ncol = facet_ncol, ...)
    return(p)
}


############################################################################
############################################################################
### Calculate the cell density
### https://www.nature.com/articles/s41467-023-40398-4
Cell_density_fun <- function(srt= NULL, meta_col= NULL, output_path= NULL, width= 2000, height= 2000, helpme= FALSE){

    # Check if helpme parameter is TRUE
    if (helpme) {
        cat("This function shows the cell density into the UMAP.\n")
        cat("Required input parameters:\n")
        cat("- srt: seurat object.\n")
        cat("- meta_col: metadata column names in seurat object (Genotype column).\n")    
        cat("- output_path: Full path for figure output.\n")            
        cat("- width: Figure width (eg. 2000).\n")
        cat("- height: Figure height (eg. 2000).\n")
        cat("- eg. Cell_density_fun(srt= srt, meta_col= meta_col, output_path= output_path, width= 5000, height= 5000, helpme= FALSE)")
        return(invisible())
    }

    # List of required packages
    required_packages <- c("LSD", "gridExtra")
    pck_check(required_packages)

    ##
    Celldencity_plot <- paste0(output_path, 'Celldensity_plot/')
    dir.create(Celldencity_plot)
    for(meta_col_val in unique(srt@meta.data[[meta_col]])){
        xmax <- max(srt@reductions[[reductionName]]@cell.embeddings[, 1])
        xmin <- min(srt@reductions[[reductionName]]@cell.embeddings[, 1])
        ymax <- max(srt@reductions[[reductionName]]@cell.embeddings[, 2])
        ymin <- min(srt@reductions[[reductionName]]@cell.embeddings[, 2])

        umap1 <- srt@reductions[[reductionName]]@cell.embeddings[, 1]
        umap2 <- srt@reductions[[reductionName]]@cell.embeddings[, 2]
        
        sel <- which(srt@meta.data[[meta_col]] == meta_col_val)
        png (paste0(Celldencity_plot,'Cell_dencity_plot_at_',meta_col,'_level_',meta_col_val,'.png'), width = width, height = height, pointsize=10, res = 300)
        heatscatter(umap1[sel], umap2[sel], add.contour = T, greyscale = F, 
                        xlab = "UMAP_1", ylab = "UMAP_2", xlim = c(xmin, xmax), 
                        ylim = c(ymin, ymax), cex = 0.3, colpal = "bl2rd", main = meta_col_val)
        dev.off()
    }
}












############################################################################
############################################################################
### Valcano plot between DE genes
Valcano_plot<- function(srt= NULL, meta_col= NULL, output_path= NULL, width= 2000, height= 2000, helpme= FALSE){


    # Check if helpme parameter is TRUE
    if (helpme) {
        cat("This function Generate the Valcano plot between the DE genes.\n")
        cat("Required input parameters:\n")
        cat("- srt: seurat object.\n")
        cat("- meta_col: metadata column names in seurat object (Genotype column).\n")    
        cat("- output_path: Full path for figure output.\n")            
        cat("- width: Figure width (eg. 2000).\n")
        cat("- height: Figure height (eg. 2000).\n")
        cat("- eg. Cell_density_fun(srt= srt, meta_col= meta_col, output_path= output_path, width= 5000, height= 5000, helpme= FALSE)")
        return(invisible())
    }

    Idents(srt) <- meta_col
    # Run differential expression analysis using FindAllMarkers
    srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    # Save the marker file
    saveRDS(srt.markers, file = paste0(output_path, "srt_top_markers_between_",meta_col,".rds"))
    #srtptwo <-readRDS(file = paste0(projDir, "srt_top_markers_between_genotype.rds"))
    srt.markers%>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC) -> top30
    write.csv(top30, paste0(output_path, "srt_top_markers_between_",meta_col,"_top30_gene.csv"), row.names=FALSE)

    # Define threshold for highly expressed genes
    logfc_threshold <- 0.25  # Adjust as desired
    pval_threshold <- 0.05  # Adjust as desired

    # Filter genes based on log fold change and p-value
    highly_expressed_genes <- srt.markers[abs(srt.markers$avg_log2FC) >= logfc_threshold & srt.markers$p_val_adj <= pval_threshold, ]

    # Define the genotypes for comparison
    genotypes <- unique(srt$Genotype_DM)

    # Create a pairwise comparison of genotypes
    genotype_pairs <- combn(genotypes, 2, simplify = FALSE)

    for (pair in genotype_pairs) {
        genotype1 <- pair[1]
        genotype2 <- pair[2]
        comparison_data <- highly_expressed_genes[highly_expressed_genes$cluster %in% c(genotype1, genotype2), ]
        
        comparison_data%>%
        group_by(cluster) %>%
        top_n(n = 20, wt = avg_log2FC) -> top_genes

        # Create the volcano plot
        volcano_plot <- ggplot(comparison_data, aes(x = ifelse(cluster == genotype1, -avg_log2FC, avg_log2FC), y = -log10(p_val_adj))) +
        geom_point(aes(color = cluster)) +
        #geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "red") +
        geom_text_repel(data = top_genes, aes(label = gene), hjust = 0, vjust = 0) +
        labs(x = "Average Fold Change", y = "-log10(Adjusted p-value)", title = paste0(genotype1, " vs. ", genotype2)) +
        theme_minimal()

        # Print or save the volcano plot
        png (paste0(output_path, 'Plots/Valcano_plot_between_',genotype1, '_vs_', genotype2,'_top_20_genes_dispplay.png'), width = width, height = height, res = 300)
        print (volcano_plot)
        dev.off()
    }

}