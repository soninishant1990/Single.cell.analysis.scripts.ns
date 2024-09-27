
######################################
library(DirichletReg)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(gridExtra)

# metaGroupName1 = 'orig.ident' ## provide sample column names
# metaGroupName2 = 'groupID' ## provide genotype or sampling column name
# metaGroupName3 = 'celltype2' ### Provide the cell type or any other metadata column name for which you want to create a boxplot. 
# ## here we can provide any name. so basically we are giving one higher name to all subtype name for metaGroupName3
# srt$immnune <- 'immnune'
# metaGroupName4 <- 'immnune'

# celltypes_pal_I1 = color.list[[metaGroupName2]] # provide list of color matching with number of metaGroupName3 column

# # FIGURE S1C Composition plots by sampling ####
# metaGroupNames = c(metaGroupName1,metaGroupName2, metaGroupName3, metaGroupName4)
# metalength <- length(metaGroupNames)

# ## Calculate the dirichlet pvalue 
# ccc_bar2 = cellComp.for.dirichlet.function1 (
#   seurat_obj = srt, 
#   metaGroups = metaGroupNames,
#   prop = TRUE,
#   ptable_factor = c(1,metalength),
#   returnDF = T
#   )


# # Generate the boxplot
# ccc_bar5 = cellComp.for.dirichlet.function2 (
#   seurat_obj = srt,
#   metaGroups = c(metaGroupName1, metaGroupName2, metaGroupName2, metaGroupName3),
#   plot_as = 'box',
#   pal = celltypes_pal_I1,
#   facet_ncol = 4,
#   pair_com = NULL, # xy.list if applicable
#   Pvalue_cal = FALSE,
#   Pvalue_method = 't_test', # "wilcox_test", "t_test", etc.
#   hide_ns = TRUE,
#   Pvalue_label = "{p.signif}", # "p", "p.adj", "p.signif", etc.
#   p.adjust.method = 'none', # "holm", "hochberg", etc.
#     Dirichlet.df = drc_res_df1
# )


# # Example with gridExtra
# library(gridExtra)

# # Save as PDF
# pdf(paste0(projDir, 'Plots/cell_composition_', metaGroupName3, '_Group_id_level_barboxplots_dirichlet.pdf'), width = 12, height = 7)

# # Arrange plots in a grid (adjust ncol as per your preference)
# grid.arrange(grobs = ccc_bar5, ncol = 4)

# # Save and close the PDF file
# dev.off()

# # Save as PNG
# png(paste0(projDir, 'Plots/cell_composition_', metaGroupName3, '_Group_id_level_barboxplots_dirichlet.png'), width = 4000, height = 1800, res = 300)

# # Arrange plots in a grid (adjust ncol as per your preference)
# grid.arrange(grobs = ccc_bar5, ncol = 4)

# # Save and close the PNG file
# dev.off()


## function for Calculating the dirichlet pvalue 
cellComp.for.dirichlet.function1 = function (
    seurat_obj = NULL, 
    metaGroups = NULL, # vector of at least 3 metaGroups e.g. c('orig.ident','celltypes','celltypes'),
    #plot_as = 'box', # box or bar 
    #pal = NULL,
    prop = TRUE,
    ptable_factor = 1, # specify which column of the data.frame or seurat object metadata should be used to compute proportions
    #facet_ncol = 20,
    #facet_scales = 'free',
    #subset_prop = NULL, # subset prop table by any group in any column
    removeNA = TRUE,
    returnDF = FALSE
    )
    {
    
    if (is.data.frame (seurat_obj))
    {
    meta_groups_df = seurat_obj[,metaGroups]
    } else {
    meta_groups_df = seurat_obj@meta.data[,metaGroups]
    }
    # Refactor to remove 0 groups
    #meta_groups_df =  as.data.frame(lapply(unclass(meta_groups_df),as.character),stringsAsFactors=T)
    #if(is.null(pal)) pal = rainbow (length(unique(meta_groups_df[,2])))
    #if(is.null(pal) & plot_as == 'box') pal = rainbow (length(unique(meta_groups_df[,3])))
    if (prop)
        {
        ccomp_df = as.data.frame (prop.table (table (meta_groups_df),ptable_factor))
        ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
        } else {
        ccomp_df = as.data.frame (table (meta_groups_df))	
        }
    if(removeNA) ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
    # if (!is.null (subset_prop)) 
    #     {
    #     subset_col = unlist(sapply (seq(ncol(ccomp_df)), function(x) if(any(ccomp_df[,x] %in% subset_prop)) colnames(ccomp_df)[x]))
    #     ccomp_df = ccomp_df[ccomp_df[,subset_col] %in% subset_prop,]
    #     }  


        ### Run Diricthlet regression and extract pvalues ####
    drc_res = list()
    i <- levels (ccomp_df[,metaGroups[[4]]])

    # Define the genotypes for comparison
    genotypes <- unique(seurat_obj@meta.data[,metaGroups[[2]]])

    # Create a pairwise comparison of genotypes
    pair_list <- combn(genotypes, 2, simplify = FALSE)

    ## Calculate the Dirichlet p value in pairwise comparision and save into drc_res list
    for (pair_list1 in pair_list){
        pair_list2 <- c(pair_list1[[1]], pair_list1[[2]])
        pair_list.names <- paste0(pair_list1[[1]],'.',pair_list1[[2]])
        ccc_bar3 <- ccomp_df[ccomp_df[[metaGroups[[2]]]] %in% pair_list2, ]
        cc_box1 = spread(ccc_bar3[ccc_bar3[[metaGroups[[4]]]] == i,], key = metaGroups[[3]], value = Freq)
        tmp = cc_box1[,4:ncol(cc_box1)]
        tmp[is.na(tmp)] = 0
        AL = DR_data(tmp)
        res = DirichReg (AL ~ as.factor (groupID), cc_box1)
        u = summary (res)

        pvals = u$coef.mat[grep('(Intercept)', rownames(u$coef.mat), invert=T), 4]
        v = names(pvals)
        #pvals1 = u$coef.mat[grep('(Intercept)', rownames(u$coef.mat), invert=T), 4]
        #v = names(pvals1)
        pvals = matrix(pvals, ncol=length(u$varnames))
        rownames(pvals) = gsub(metaGroups[[2]], '', v[1:nrow(pvals)])
        colnames(pvals) = u$varnames
        pvals = as.data.frame (pvals)
        pvals$groupID = rownames (pvals)
        pvals = gather (pvals, celltype, p, 1:(ncol(pvals)-1))
        pvals$groupID = gsub ('as.factor\\(\\)','',pvals$groupID)
        pvals$group1 = pair_list1[[1]]
        pvals[, metaGroups[[2]]] = NULL
        pvals$group2 = pair_list1[[2]]
        pvals$p_adj = p.adjust (pvals$p, method = 'fdr')
        pvals$p_sig = ''
        pvals$p_sig = ifelse (pvals$p_adj <= 0.05, ifelse (pvals$p_adj <= 0.01,ifelse (pvals$p_adj <= 0.001,'***','**'),'*'),'ns')
        df <- pvals
        drc_res[[pair_list.names]] = df
    }

    ### Merge dataframe list inot one
    drc_res_df = do.call (rbind, drc_res)
    drc_res_df1 <- drc_res_df[order(drc_res_df$celltype), ]

    ### Calculate the y axis position in plot to show asterisk or pvalue at correct place
    ccc_bar4 = split (ccomp_df, ccomp_df[,metaGroups[[3]]])

    y_max = do.call(rbind, lapply(ccc_bar4, function(x) {
      tpm = boxplot(x$Freq ~ x[[metaGroups[[2]]]])$stats
      tpm = max(tpm[5,])
      tmp = data.frame(
        celltype = x[[metaGroups[[2]]]],
        groupID = x[[metaGroups[[3]]]],
        y.position = tpm
      )
      tmp = tmp[!duplicated(tmp$celltype), ]
      tmp
    }))
    # Ensure `rownames` are reset
    rownames(y_max) <- NULL

    ## Add y position into main dataframe
    drc_res_df1$y.position = y_max$y.position 
    ## Add pairwise value in immnune column
    drc_res_df1$immune = gsub ('\\.\\d','',rownames (drc_res_df1))
    # remove non significant value from datafram
    drc_res_df1 <- subset(drc_res_df1, p_sig != 'ns')

    return(drc_res_df1)
}




# Function for Generate the boxplot
cellComp.for.dirichlet.function2  = function (
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
                    p.adjust.method = NULL,
                    Dirichlet.df = NULL
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
       ### generated plot for each celltype annd add pvalue      
                                
        boxplot.list = list()
      for (celltypename in unique(ccomp_df[,metaGroups[length(metaGroups)]])){
        ccomp_df1 <- ccomp_df[ccomp_df[, metaGroups[length(metaGroups)]] == celltypename, ]
        colnames (ccomp_df1) = c(paste0('Var_',seq_along(metaGroups)), 'proportion')
        box_p = ggplot (ccomp_df1, aes (x= Var_2, y= proportion)) +
          geom_boxplot(aes (fill= Var_3)) +
          theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          scale_fill_manual (values= pal) +
          ggtitle(celltypename) +
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
                strip.background = element_rect(fill = "transparent")) # Make facet strip background transparent
    #box_p
        ## add pvalue
          #boxplot.list <- list[]
        #drc_res_df1_1 <- subset(drc_res_df1, celltype == celltypename)
        drc_res_df1_1 <-   Dirichlet.df[Dirichlet.df[, 'celltype'] == celltypename, ]
        # Add p-values using stat_pvalue_manual
        ccc_bar6 = box_p + stat_pvalue_manual(
          drc_res_df1_1,
          step.increase = 0.05,
          y.position = 'y.position',
          #x = c("group1", "group2"), # Adjusted to match your grouping variable
          hide.ns = TRUE,
          color = 'grey20',
          label = "p_sig",
          vjust = 0.5,
          tip.length = 0.02
        )
        boxplot.list[[celltypename]] <- ccc_bar6
        #ccc_bar6
        #print(celltypename)
    }

  return (boxplot.list)
}