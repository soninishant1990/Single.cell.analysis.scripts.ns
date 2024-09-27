
### Usefull funtion for scRNA analysis 

############################################################################
############################################################################
## Cell composition function
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
    geom_boxplot(aes (fill= Var_3), outlier.shape = NA) +
     geom_jitter (color="black", size=0.8, alpha=0.9) +
    theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual (values= pal) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),  # Set x-axis text to bold
        axis.text.y = element_text(face = "bold"),  # Set y-axis text to bold
        axis.title = element_text(face = "bold"),  # Set axis titles to bold
        strip.text = element_text(face = "bold"),  # Set facet strip text to bold
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.background = element_blank(),  # Remove legend background
        legend.box.background = element_blank(),  # Remove legend box background
        legend.key = element_blank(),  # Remove legend keys (color squares)
        legend.title = element_blank(),  # Remove legend title background
        strip.background = element_rect(fill = "transparent")) #+  # Make facet strip background transparent
  #guides(fill = guide_legend(override.aes = list(fill = NA)))
    
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
    #theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual (values= pal) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),  # Set x-axis text to bold
        axis.text.y = element_text(face = "bold"),  # Set y-axis text to bold
        axis.title = element_text(face = "bold"),  # Set axis titles to bold
        strip.text = element_text(face = "bold"),  # Set facet strip text to bold
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.line = element_line(color = "black")
         )
    if (length(metaGroups) == 3) box_p = box_p + facet_wrap (~Var_3, scales=facet_scales, ncol=facet_ncol)
    }
return (box_p)
}


############################################################################
############################################################################
### TCGA data analysis function

### Funtion to calculate expression for each modules
Add_celltype_expression <- function (x, gene_list) {

    # Find the maximum length of the elements in gene_list
    max_length <- max(sapply(gene_list, length))

    # Pad the shorter elements with NA values
    for (i in 1:length(gene_list)) {
      if (length(gene_list[[i]]) < max_length) {
        gene_list[[i]] <- c(gene_list[[i]], rep(NA, max_length - length(gene_list[[i]])))
      }
    }

    # Convert each element to a data frame
    for (i in 1:length(gene_list)) {
      gene_list[[i]] <- as.data.frame(t(as.data.frame(gene_list[[i]])))
    }

    # Create a data frame with row names
    genes.df <- as.data.frame(do.call(rbind, gene_list))
    genes.df$GeneSet <- names(gene_list)
    
    
    
    genes.df <- rbindlist(gene_list)
    genes.df <- as.data.frame(genes.df)
    rownames(genes.df) <- names(gene_list)
    gmt.program.list <- list()

    for (i in 1:nrow(genes.df)) {
      row <- genes.df[i,]
      row <- as.data.frame(t(as.data.frame(row)))
      row$cluster <- i-1
      row$Annotation <- rownames(genes.df)[i]
      colnames(row)[1] <- "gene"
      gmt.program.list[[i]] <- row
    }

    modules_df <- as.data.frame(rbindlist(gmt.program.list, fill = T))
    modules_df <- modules_df[complete.cases(modules_df), ]
    modules <- names(gene_list)

    top10genes <- modules_df%>% group_by(cluster)
    top.genes<-list()

    for (i in 1:length(modules)){
      name<-modules[i]
      genes<-as.character(top10genes$gene[which(top10genes$cluster %in% (i-1))]) # for celltypeIM, etc
      top.genes[[name]]<-genes
    }

    modules <- names(top.genes)

    # Some more reformatting
    input.df <- as.matrix(log10(GetAssayData(gbm.2018,slot='counts')+1))
    rownames(input.df)<-rownames(GetAssayData(gbm.2018,slot='counts'))

    exprs.ave.list<-list()

    for (i in 1:length(modules)){
      markers<-toupper(top.genes[[i]])
      module_name<-modules[i]
      print(module_name)
      index<-which(rownames(input.df) %in% markers)
      print(length(index))
      exprs<-input.df[index,, drop = F]
      rownames(exprs)<-rownames(input.df)[index]

      if (length(index) == 1){
        exprs.ave.list[[i]]<-exprs
      }
      else if (length(index) > 1){
        exprs.ave.list[[i]]<-colMeans(exprs)
      }
      else {exprs.ave.list[[i]]<-NA
      print(module_name)}
    }

    modules<-modules[which(!is.na(exprs.ave.list))]
    top.genes<-top.genes[which(!is.na(exprs.ave.list))]
    exprs.ave.list[which(is.na(exprs.ave.list))]<-NULL
    exprs.ave.df<-as.data.frame(t(do.call(rbind,exprs.ave.list)))
    names(exprs.ave.df) <- modules

    tmp<-cbind(gbm.2018@meta.data, exprs.ave.df)

    # Looking at the number of each subtype
    #table(tmp$new_subtype_basedon_mut)

    #dplyr::count(tmp, new_subtype_basedon_mut, sort = TRUE)
    return(list(tmp, modules))
}

### Function to get plot grid 
find_grid_side_length <- function(number) {
  sqrt_val <- sqrt(number)
  floor_sqrt <- floor(sqrt_val)
  
  # Check if the number is a perfect square
  if (sqrt_val == floor_sqrt) {
    return(floor_sqrt)
  } else {
    # If not a perfect square, return the nearest larger square's side length
    return(floor_sqrt + 1)
  }
}
### Boxplot function
boxplot.fun1 <- function(obj = NULL, col.name = NULL, genotype.color = NULL, level.factor = NULL, file.name = NULL,
                       width = NULL, height = NULL){
    obj <- subset(obj, get(col.name) != 'na')
    obj[,col.name] <- factor(obj[,col.name], levels = level.factor)
    
    ### pair for p value comparision
    pair_list <- lapply(2:length(level.factor), function(k) combn(level.factor, k, simplify = FALSE))
    pair_list <- unlist(pair_list, recursive = FALSE)
    # Filter out combinations with more than two elements
    pair_list <- pair_list[sapply(pair_list, function(x) length(x) == 2)]
    my_comparisons <- pair_list
                                  
    tmp_table <- table(obj[,col.name])
    tmp_list <- as.list(tmp_table)
    subtype_values <- paste(names(tmp_list), tmp_list, sep = " = ")
    subtype_values_combined <- paste(subtype_values, collapse = ", ")
    main_title <- paste0("These subtype specify through specific mutation ",subtype_values_combined, ". Only mut patient selected from Brennan paper assign subtype.")                       
                                  
    mean.plots <- list()                              
    for (module in modules){
      print(module)

      mean.plot <- ggboxplot(obj, x = col.name, y = module,
                             color = col.name, palette = genotype.color,
                             title = module) +
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = TRUE) +
      rotate_x_text(90) +
      theme(plot.title = element_text(size = 8, face = "bold"), legend.position = "none") +
         labs(x = NULL)  # Remove y-axis title

      mean.plots[[module]] <- mean.plot
    }
    ## get plot grid
    number <- length(modules)
    grid_side_length <- find_grid_side_length(number)
    # Combine all plots into a single plot grid
    combined_plot <- cowplot::plot_grid(plotlist = mean.plots, ncol = grid_side_length)

    # Add main text at the top of the plot
    combined_plot <- combined_plot +
      ggtitle(main_title) +
      theme(plot.title = element_text(size = 7, face = "bold", margin = margin(b = 20)))
                                  
    # Save as PDF
    pdf(file.name, width = width, height = height, useDingbats=FALSE)
    print(combined_plot)
    dev.off()
                                  
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
### cellphonedb analysis
message('Cellphonedb_ds_analysis function to generate ligand receptor dataframe')
message('Cellphonedb_dotplot function to generate ligand receptor dotplot by pairwise comparsion between genotype')

Cellphonedb_ds_analysis <- function(input.dir=NULL, out.dir = NULL, sample.names = NULL, sample_col_name = NULL, genotype_colname =NULL,
                                   genotype_analysis = NULL){
    ### Read all singnificant file and add into single datafram list
    results.df.list <- list()
    for (i in 1:length(unique(sub_srt@meta.data[,sample_id_col]))) {
      signif.means.df <- read_delim(paste0(input.dir, "/", unique(sub_srt@meta.data[,sample_id_col])[i], "/significant_means.txt"), col_names =TRUE, show_col_types = FALSE, delim = '\t')
      # Get the current column names
      current_names <- colnames(signif.means.df)
      # Replace dots with underscores in column names
      new_names <- gsub("\\.", "_", current_names)
      # Assign the modified column names back to the dataframe
      colnames(signif.means.df) <- new_names
      info_columns <- signif.means.df[,c(1,2,5,6,12)]
      cell_pair_columns <- signif.means.df[,13:ncol(signif.means.df)]

      signif.means.df <- cbind(info_columns, cell_pair_columns)
      results.df.list[[i]] <- signif.means.df
    }

    names(results.df.list) <- unique(sub_srt@meta.data[,sample_id_col])
    results.df.list <- results.df.list[unique(sub_srt@meta.data[,sample_id_col])]
    
    # find the union of all interactions in your dataset
    combined_df <- rbindlist(results.df.list, fill=TRUE)
    combined_df <- combined_df[!duplicated(combined_df$id_cp_interaction)]
    combined_df <- as.data.frame(combined_df)
    rownames(combined_df) <- combined_df$id_cp_interaction
    
    # this should make all the data frames have the same columns
    full_columns <- setNames(data.frame(matrix(ncol = ncol(combined_df), nrow = 0)), colnames(combined_df)) # has all column names

    for (i in 1:length(results.df.list)) {
      sample <- results.df.list[[i]]
      sample <- as.data.frame(rbindlist(list(full_columns, sample), fill = T))
      sample <- as.data.frame(sample)
      sample <- sample[!duplicated(sample$id_cp_interaction),]
      rownames(sample) <- sample$id_cp_interaction
      results.df.list[[i]] <- sample
    }

    # this should make all the data frames have the same rows
    full_rows <- setNames(data.frame(matrix(ncol = nrow(combined_df), nrow = 0)), combined_df$id_cp_interaction) # has all row names

    for (i in 1:length(results.df.list)) {
      sample <- as.data.frame(results.df.list[[i]])
      # save cell pair names
      cols <- colnames(sample)

      sample <- as.data.frame(t(sample)) # temporarily transpose matrix to rbind full row names
      sample <- as.data.frame(rbindlist(list(full_rows, sample), fill = T))

      sample <- as.data.frame(t(sample))
      sample[,1:4] <- combined_df[,1:4] # fills in empty rows due to extending the data frame
      # reassign cell pair names
      colnames(sample) <- cols
      results.df.list[[i]] <- sample
    }

    mut_signif <- data.frame(matrix(0, nrow = nrow(combined_df), ncol = ncol(combined_df), dimnames = list(rownames(combined_df), colnames(combined_df))))
    wt_signif <- data.frame(matrix(0, nrow = nrow(combined_df), ncol = ncol(combined_df), dimnames = list(rownames(combined_df), colnames(combined_df))))
    diff_prop <- data.frame(matrix(0, nrow = nrow(combined_df), ncol = ncol(combined_df), dimnames = list(rownames(combined_df), colnames(combined_df))))
    fishertest_pvals <- data.frame(matrix(0, nrow = nrow(combined_df), ncol = ncol(combined_df), dimnames = list(rownames(combined_df), colnames(combined_df))))
    
    ## Get list of sample sepratly for both condition
    qun_srt = sub_srt@meta.data[!duplicated(sub_srt@meta.data[,sample_col_name]),]
    qun_srt1 <- qun_srt[,c(sample_col_name,genotype_colname)]
    qun_srt_type_one <- qun_srt1[qun_srt1[genotype_colname] == genotype_analysis[1],]
    qun_srt_type_one_sample_list <- qun_srt_type_one[,sample_col_name]
    qun_srt_type_two <- qun_srt1[qun_srt1[genotype_colname] == genotype_analysis[2],]
    qun_srt_type_two_sample_list <- qun_srt_type_two[,sample_col_name]
    
    
    ### Here type one (qun_srt_type_one_sample_list) is wt and type two (qun_srt_type_two_sample_list) is mut
    ### 
    for (i in 1:length(qun_srt_type_two_sample_list)) {
      sample <- results.df.list[[qun_srt_type_two_sample_list[i]]]
      for (j in 6:ncol(sample)) {
        for (k in 1:nrow(sample)) {
          value <- sample[k,j]
          if (!is.na(value)) {
            wt_signif[k,j] <- wt_signif[k,j] + 1
          }
        }
      }
    }

    ### 
    for (i in 1:length(qun_srt_type_one_sample_list)) {
      sample <- results.df.list[[qun_srt_type_one_sample_list[i]]]
      for (j in 6:ncol(sample)) {
        for (k in 1:nrow(sample)) {
          value <- sample[k,j]
          if (!is.na(value)) {
            mut_signif[k,j] <- mut_signif[k,j] + 1
          }
        }
      }
    }

    ###
    for (j in 6:ncol(diff_prop)) {
      for (k in 1:nrow(diff_prop)) {
        diff_prop[k,j] <- mut_signif[k,j]/length(qun_srt_type_one_sample_list) - wt_signif[k,j]/length(qun_srt_type_two_sample_list)
      }
    }

    ###
    for (j in 6:ncol(fishertest_pvals)) {
      for (k in 1:nrow(fishertest_pvals)) {
        dat <- data.frame(
          "signif_no" = c(length(qun_srt_type_one_sample_list)-mut_signif[k,j], length(qun_srt_type_two_sample_list)-wt_signif[k,j]),
          "signif_yes" = c(mut_signif[k,j], wt_signif[k,j]),
          row.names = genotype_analysis,
          stringsAsFactors = FALSE
        )
        test <- fisher.test(dat)
        fishertest_pvals[k,j] <- -log10(test$p.value)
      }
    }
    
    mut_signif_vector <- as.data.frame(t(as.data.frame(unmatrix(as.matrix(mut_signif[,6:ncol(mut_signif)])))))
    wt_signif_vector <- as.data.frame(t(as.data.frame(unmatrix(as.matrix(wt_signif[,6:ncol(wt_signif)])))))
    diff_prop_vector <- as.data.frame(t(as.data.frame(unmatrix(as.matrix(diff_prop[,6:ncol(diff_prop)])))))
    fishertest_pvals_vector <- as.data.frame(t(as.data.frame(unmatrix(as.matrix(fishertest_pvals[,6:ncol(fishertest_pvals)])))))

    df <- as.data.frame(t(rbindlist(list(mut_signif_vector, wt_signif_vector, diff_prop_vector, fishertest_pvals_vector))))
    colnames(df) <- c("mut_num_signif", "wt_num_signif", "mut_wt_prop_diff", "mut_wt_neglogpval_fishertest")
    df$int_name <- sapply(strsplit(rownames(df),":"), `[`, 1)
    df$cell.pair <- sapply(strsplit(rownames(df),":"), `[`, 2)
    df$cell.a <- sapply(strsplit(df$cell.pair,"\\."), `[`, 1)
    df$cell.b <- sapply(strsplit(df$cell.pair,"\\."), `[`, 2)
    
    
    interaction_input_path <- "/ahg/regevdata/projects/ICA_Lung/Nishant/Colon_cancer_atlas/5_Cellphonedb/cpdb/"
    load(file = paste0(interaction_input_path, "interaction_table_updated.rda"))
    interaction <- interaction_table2
    interaction <- as.data.frame(interaction)
    
    sum(duplicated(interaction$id_cp_interaction))
    rownames(interaction) <- interaction$id_cp_interaction
    interaction2 <- interaction[df$int_name,]
    interaction2$gene_a <- interaction2$multidata_1_id
    interaction2$gene_b <- interaction2$multidata_2_id
    interaction2 <- interaction2[,c("id_cp_interaction", "gene_a", "gene_b")]
    lig_rec_pair_names <- combined_df[,c(1:2)]
    lig_rec_pair_names2 <- lig_rec_pair_names[df$int_name,]
    all(rownames(lig_rec_pair_names2) == rownames(interaction2))
    interaction2$interacting_pair <- lig_rec_pair_names2$interacting_pair
    interaction2 <- interaction2[2:ncol(interaction2)]
    df2 <- cbind(df, interaction2)
    
    lig_rec_differential_highlvl_df <- df2 #df2
    save(lig_rec_differential_highlvl_df, file = paste0(out.dir , "/lig_rec_differential_highlvl_df.Rda"))
    write.csv(as.matrix(lig_rec_differential_highlvl_df), file = paste0(out.dir , "/lig_rec_differential_highlvl_df.csv"))
    range(lig_rec_differential_highlvl_df$mut_wt_neglogpval_fishertest, na.rm = T)
    sum(lig_rec_differential_highlvl_df$mut_wt_neglogpval_fishertest == 0, na.rm = T)
    sum(is.na(lig_rec_differential_highlvl_df$mut_wt_neglogpval_fishertest))
    sum(lig_rec_differential_highlvl_df$mut_num_signif == 0 & lig_rec_differential_highlvl_df$wt_num_signif == 0, na.rm = T)
    lig_rec_differential_highlvl_df_not_allzeros <- lig_rec_differential_highlvl_df[!(lig_rec_differential_highlvl_df$mut_num_signif == 0 & lig_rec_differential_highlvl_df$wt_num_signif == 0),]
    save(lig_rec_differential_highlvl_df_not_allzeros, file = paste0(out.dir, "/lig_rec_differential_highlvl_df.Filtered.Rda"))
    write.csv(as.matrix(lig_rec_differential_highlvl_df_not_allzeros), file = paste0(out.dir, "/lig_rec_differential_highlvl_df.Filtered.csv"))
    load(file = paste0(out.dir , "/lig_rec_differential_highlvl_df.Filtered.Rda"))
    
}




############################################################################
############################################################################
### Neftel scatter plot

# Generate GBM cell states as in paper Nefter et al
# On y axis = D = max(SCopc,SCnpc) - max(SCac,SCmes)
ModScoreCor = function (seurat_obj, geneset_list, listName, cor_threshold = NULL, pos_threshold = .1, outdir)
        {        
        message ('Run AddModuleScore')
        seurat_obj = AddModuleScore (seurat_obj, geneset_list)
        seurat_obj@meta.data = seurat_obj@meta.data[, !colnames (seurat_obj@meta.data) %in% names (geneset_list)]
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %in% paste0('Cluster',seq_along(geneset_list))] = names (geneset_list)
        message (paste('Annotate cells based on highest module score and store in column:',paste0(listName, '_r',cor_threshold,'_max')))
        if (length (geneset_list) == 1) 
            {
            seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = ifelse (seurat_obj@meta.data[,names (geneset_list)] > pos_threshold, 'pos','neg')
            pdf (paste0(outdir, listName, '_modulescore_distribution_cor_threshold_',cor_threshold,'_score_',pos_threshold,'.pdf'))
            hist (seurat_obj@meta.data[,names(geneset_list)])
          abline (v = pos_threshold)
          dev.off()
            } else {
            seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = sapply (seq_along(colnames(seurat_obj)), function(x) colnames(seurat_obj@meta.data[,names(geneset_list)])[which.max (seurat_obj@meta.data[x,names(geneset_list)])])               
            }
        if (!is.null(cor_threshold))
                {
                message ('cor_threshold provided! Filtering gene sets based on initial correlation to module score')    
                filtered_geneset_list = list()
                geneset_cor_list = list()
                for (i in names(geneset_list))
                        {       
                        geneset_cor = cor (seurat_obj@meta.data[,i], as.matrix(t(seurat_obj@assays$RNA@data[rownames(seurat_obj@assays$RNA@data) %in% geneset_list[[i]],])))
                        geneset_cor_list[[i]] = geneset_cor
                        geneset_cor_names = colnames (geneset_cor)[geneset_cor > cor_threshold]
                        geneset_cor_names = geneset_cor_names[!is.na (geneset_cor_names)]
                        filtered_geneset_list[[i]] = geneset_cor_names
                        }
                if (!is.null (outdir)) 
                        {
                        lapply (seq_along(filtered_geneset_list), function(x) write.csv (filtered_geneset_list[[x]], paste0(outdir,'Corfiltered_Module_score_gene_list_', names(filtered_geneset_list)[x],'.csv')))
                        pdf (paste0(outdir, listName, 'Corfiltered_modulescore_distribution.pdf'))
                        lapply (seq_along(filtered_geneset_list), function(x) 
                                {
                                hist (geneset_cor_list[[x]], title = names(geneset_cor_list)[x])
                                abline (v = cor_threshold)
                                })
                        dev.off()
                        }
                message ('Re-run AddModuleScore using corfiltered genes')
                seurat_obj = AddModuleScore (seurat_obj, filtered_geneset_list, name = listName)
                seurat_obj@meta.data = seurat_obj@meta.data[, !colnames (seurat_obj@meta.data) %in% paste0(names(geneset_list),'_r',cor_threshold)]
                colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %in% paste0('Cluster',seq_along(geneset_list))] = paste0(names(geneset_list),'_r',cor_threshold)
                if (length (geneset_list) == 1) 
                    {
                            seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = ifelse (seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)] > pos_threshold, 'pos','neg')
                            pdf (paste0(outdir, listName, '_modulescore_distribution_cor_threshold_',cor_threshold,'_score_',pos_threshold,'.pdf'))
                            hist (seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)])
                        abline (v = pos_threshold)
                        dev.off()
                            } else {
                    seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = sapply (seq_along(colnames(seurat_obj)), function(x) colnames(seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)])[which.max (seurat_obj@meta.data[x,paste0(names(geneset_list),'_r',cor_threshold)])])        
                    }
                }
        return (seurat_obj)
        } 


Nefstates = function (seurat_obj = srt, mouse=T)
    {
    require (Seurat)    
    message ('Read csv and compute module scores')
    if (mouse) gbm_modules_6states = readRDS (paste0 ('../../Genotype.prj/data/mouse_GBM_netfel_6_states.rds')) else
    gbm_modules_6states = readRDS (paste0 ('../../Genotype.prj/data/human_GBM_netfel_6_states.rds'))        
    #source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/cellcycle.R')
    message ('Cell cycle analysis')
    if (mouse == TRUE) {
        m.g2m.genes = readRDS ('../../Genotype.prj/data/cellcycle_mouse.g2m.genes.rds')
        m.s.genes = readRDS ('../../Genotype.prj/data/cellcycle_mouse.s.genes.rds')
        s.genes <- m.s.genes
        g2m.genes <- m.g2m.genes
    } else if (org == 'human'){
        s.genes <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes
    }
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
    seurat_obj$cc_score = seurat_obj$G2M.Score + seurat_obj$S.Score
    
    
    seurat_obj = ModScoreCor (
        seurat_obj = seurat_obj, 
        geneset_list = gbm_modules_6states, 
        cor_threshold = NULL, 
        pos_threshold = -.03,
        listName = 'neftel6', 
        outdir = NULL)
    seurat_obj_df = as.data.frame (seurat_obj@meta.data)
    seurat_obj_df$MESlike = apply (seurat_obj_df[,c('MES1','MES2')],1,max)
    seurat_obj_df$NPClike = apply (seurat_obj_df[,c('NPC1','NPC2')],1,max)
    
    message ('compute y axis')
    max_AC_or_MES = pmax(seurat_obj_df$AC, seurat_obj_df$MESlike)
    max_OPC_or_NPC = pmax(seurat_obj_df$OPC, seurat_obj_df$NPClike)
    seurat_obj_df$y_axis <- log2(abs(max_OPC_or_NPC - max_AC_or_MES) + 1)
    seurat_obj_df$y_axis[max_AC_or_MES > max_OPC_or_NPC] <- -1 * seurat_obj_df$y_axis[max_AC_or_MES > max_OPC_or_NPC]

    message ('compute x axis')
    seurat_obj_df$x_axis = 0
    seurat_obj_df$x_axis[seurat_obj_df$y_axis > 0] = log2(abs(seurat_obj_df$OPC - seurat_obj_df$NPClike) + 1)[seurat_obj_df$y_axis > 0]
    seurat_obj_df$x_axis[seurat_obj_df$y_axis > 0 & seurat_obj_df$OPC > seurat_obj_df$NPClike] <- -1 * seurat_obj_df$x_axis[seurat_obj_df$y_axis > 0 & seurat_obj_df$OPC > seurat_obj_df$NPClike]
    
    seurat_obj_df$x_axis[seurat_obj_df$y_axis < 0] = log2(abs(seurat_obj_df$AC - seurat_obj_df$MESlike) + 1)[seurat_obj_df$y_axis < 0]
    seurat_obj_df$x_axis[seurat_obj_df$y_axis < 0 & seurat_obj_df$AC > seurat_obj_df$MESlike] <- -1 * seurat_obj_df$x_axis[seurat_obj_df$y_axis < 0 & seurat_obj_df$AC > seurat_obj_df$MESlike]
        
    message ('assign states to cells')
    seurat_obj_df$Tstate=0
    seurat_obj_df$Tstate[seurat_obj_df$x_axis < 0 & seurat_obj_df$y_axis >= 0] = 'OPC-like' 
    seurat_obj_df$Tstate[seurat_obj_df$x_axis >= 0 & seurat_obj_df$y_axis > 0] = 'NPC-like' 
    seurat_obj_df$Tstate[seurat_obj_df$x_axis < 0 & seurat_obj_df$y_axis < 0] = 'AC-like' 
    seurat_obj_df$Tstate[seurat_obj_df$x_axis > 0 & seurat_obj_df$y_axis < 0] = 'MES-like'
    
    return (seurat_obj_df)
    }