
## To generate ligand receptor interaction dataframe from cellphonedb result. 
## This code is use for doing pairwise comparision between groups
#pair_list <- list('EGFR.WT'=c('NF1.WT'), 'EGFR.WT'=c('PDGFB.WT'), 'NF1.WT'=c('PDGFB.WT'), 'PDGFB.IDHmut.recur' = c('PDGFB.WT.recur'))

# for (pair_name in names(pair_list)){
#     print(pair_name)
#     print(pair_list[[pair_name]])
    
#     new_pair_list <- c(pair_name, pair_list[[pair_name]])
    
#     n = 0
#     for (group_name in new_pair_list){
#         if (n==0){
#             sub_srt <- subset(srt, groupID == group_name)
#         } else{
#             sub_srt1 <- subset(srt, groupID == group_name)
#             sub_srt <- merge(sub_srt, sub_srt1)
#         }
#         n = n+1
#     }
#     print(table(sub_srt$groupID))
#     sub_srt$genotype_group = ifelse (sub_srt$groupID == pair_name,pair_name,'rest')
    
#     input.dir <- paste0(subDir,'out_test/')
#     sample.names <- list.files(paste0(subDir,'out_test/'))
#     sample_col_name <- 'orig.ident'
#     genotype_colname <- 'genotype_group' ## 
#     genotype_analysis <- c(pair_name,'rest') ## If
#     sub_srt <- sub_srt

#     out.dir <- paste0(out_file,'/',pair_name,'_vs_',paste(pair_list[[pair_name]], collapse = '_'),'_cpdb_ds_analysis/')
#     dir.create(out.dir)
    
#     ### Generate interaction file for generating dotplot
#     Cellphonedb_ds_analysis(input.dir=input.dir, out.dir = out.dir, sample.names = sample.names, 
#                                         sample_col_name = sample_col_name, genotype_colname =genotype_colname,
#                                         genotype_analysis = genotype_analysis)
    
# }


### To gene dotplot
# for (pair_name in names(pair_list)){
#     print(pair_name)
#     print(pair_list[[pair_name]])
    
#     new_pair_list <- c(pair_name, pair_list[[pair_name]])
    
#     out.dir <- paste0(out_file,'/',pair_name,'_vs_',paste(pair_list[[pair_name]], collapse = '_'),'_cpdb_ds_analysis/')
#     dir.create(out.dir)
    
    
#     ### Generate dotplot
#     input_file_name <- paste0(out.dir, "/lig_rec_differential_highlvl_df_not_allzeros.Rda")
#     ## 1 = combined
#     ## 2 = seperate dotplot for all celltype
#     ## 3 = compare between specific celltype
#     ## 4 = Compare bewteen Tumor cells vs other cells
#     dotplot_type_n = 1

#     # Select the minimum prop diff both side in mutant and wild type
#     mut_wt_prop_diff_max = 0.2
#     mut_wt_prop_diff_min = -0.2
#     genotype_colname <- paste0(pair_name,'_vs_',paste(pair_list[[pair_name]], collapse = '_')) ## This name will show into the plot
#     celltype_a = NULL
#     celltype_b = NULL

#     ## Filter out sample
#     #total_sample <- 9
#     ## If total sample is less that 20 that nu_significant_sample value read as total singnificant sample in each side other wise it will calculate percentage of the sample
#     #wt_nu_significant_sample = 2 ### min 5% sample should be significant both side mut and wt
#     #mut_nu_significant_sample = 4


#     Cellphonedb_dotplot(input_file_name = input_file_name, out.dir = out.dir, dotplot_type = dotplot_type_n, celltype_a = NULL, celltype_b = NULL,
#                         mut_wt_prop_diff_max=mut_wt_prop_diff_max, mut_wt_prop_diff_min =mut_wt_prop_diff_min,
#                        genotype_colname = genotype_colname, tumor_cell_name=NULL)

#     dotplot_type_n = 2
#     Cellphonedb_dotplot(input_file_name = input_file_name, out.dir = out.dir, dotplot_type = dotplot_type_n, celltype_a = NULL, celltype_b = NULL,
#                         mut_wt_prop_diff_max=mut_wt_prop_diff_max, mut_wt_prop_diff_min =mut_wt_prop_diff_min,
#                        genotype_colname = genotype_colname, tumor_cell_name=NULL)
    
#     dotplot_type_n = 4
#     tumor_cell_name = 'Tumor'
#     Cellphonedb_dotplot(input_file_name = input_file_name, out.dir = out.dir, dotplot_type = dotplot_type_n, celltype_a = NULL, celltype_b = NULL,
#                         mut_wt_prop_diff_max=mut_wt_prop_diff_max, mut_wt_prop_diff_min =mut_wt_prop_diff_min,
#                        genotype_colname = genotype_colname, tumor_cell_name=tumor_cell_name)
    
# }




library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(paletteer)
library (RColorBrewer)
library(hash)
library(tidyverse)
library(dplyr)
library(data.table)
library(rlist)
library(scales)
library(tidyverse)
library(gdata)

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
    save(lig_rec_differential_highlvl_df_not_allzeros, file = paste0(out.dir, "/lig_rec_differential_highlvl_df_not_allzeros.Rda"))
    write.csv(as.matrix(lig_rec_differential_highlvl_df_not_allzeros), file = paste0(out.dir, "/lig_rec_differential_highlvl_df_not_allzeros.csv"))
    load(file = paste0(out.dir , "/lig_rec_differential_highlvl_df_not_allzeros.Rda"))
    
}

Cellphonedb_dotplot <- function(input_file_name = NULL, out.dir = NULL, dotplot_type = NULL, celltype_a = NULL, 
                                celltype_b = NULL, genotype_colname = NULL, tumor_cell_name=NULL, 
                                celltype_list = NULL, Type_of_selcted_cells = NULL, top_ligand_receptor = NULL, prog_list=NULL){
    if (dotplot_type == 1){
        # To generate combined dotplot for all celltype
        
        df <- get(load(input_file_name))
        df$signif <- 0
        df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
        df$signif <- factor(df$signif, levels = c(0,1))
        select_best_pvalue = df[df$mut_wt_neglogpval_fishertest > 0.99,]

        if (!is.null(top_ligand_receptor)) {
          # Select the top 10 rows based on 'mut_wt_neglogpval_fishertest'
          top_10_rows <- head(select_best_pvalue[order(select_best_pvalue$mut_wt_prop_diff, decreasing = TRUE), ], top_ligand_receptor)
          # Select the bottom 10 rows based on 'mut_wt_neglogpval_fishertest'
          bottom_10_rows <- tail(select_best_pvalue[order(select_best_pvalue$mut_wt_prop_diff, decreasing = TRUE), ], top_ligand_receptor)
          # Combine the top and bottom rows into a single DataFrame
          select_best_pvalue <- rbind(top_10_rows, bottom_10_rows)
          output_filename <- paste0("/dotplot_lig_rec_dt_",dotplot_type,"_top_",top_ligand_receptor,".pdf")
      } else {
        output_filename <- paste0("/dotplot_lig_rec_dt_",dotplot_type,"_all_sig.pdf")
      }

        high_freq_neglogpval_fishertest = names (table (select_best_pvalue$interacting_pair))
        df = df[df$interacting_pair %in% unique(c(high_freq_neglogpval_fishertest)), ]
        dotplot_fuc(df, high_freq_neglogpval_fishertest, out.dir, output_filename, genotype_colname)  
        
    } else if (dotplot_type == 2){
        ## To make seprate dotplot for all celltype (each dotplot will show 1 celltype vs all other celltype)
        df <- get(load(input_file_name))
        for(i in unlist(dimnames(table(df$cell.a)))){
          print(i)
          
          df <- get(load(input_file_name))
          df <- lig_rec_differential_highlvl_df_not_allzeros[
            grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
            grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE), 
          ]
          df$signif <- 0
          df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
          df$signif <- factor(df$signif, levels = c(0,1))
          select_best_pvalue = df[df$mut_wt_neglogpval_fishertest > 0.99,]

          if (!is.null(top_ligand_receptor)) {
              # Select the top 10 rows based on 'mut_wt_neglogpval_fishertest'
              top_10_rows <- head(select_best_pvalue[order(select_best_pvalue$mut_wt_prop_diff, decreasing = TRUE), ], top_ligand_receptor)
              # Select the bottom 10 rows based on 'mut_wt_neglogpval_fishertest'
              bottom_10_rows <- tail(select_best_pvalue[order(select_best_pvalue$mut_wt_prop_diff, decreasing = TRUE), ], top_ligand_receptor)
              # Combine the top and bottom rows into a single DataFrame
              select_best_pvalue <- rbind(top_10_rows, bottom_10_rows)
              output_filename <- paste0("/dotplot_lig_rec_dt_",dotplot_type,"_top_",top_ligand_receptor,"_",i,".pdf")
          } else {
            output_filename <- paste0("/dotplot_lig_rec_dt_",dotplot_type,"_all_sig_",i,".pdf")
          }

          high_freq_neglogpval_fishertest = names (table (select_best_pvalue$interacting_pair))
          df = df[df$interacting_pair %in% unique(c(high_freq_neglogpval_fishertest)), ]
          
          # Get the index that sorts the 'cell.a' column based on the presence of 'CAFs'
          index <- order(grepl(i, df$cell.a, fixed = TRUE), decreasing = TRUE)
          # Reorder the DataFrame based on the index
          df <- df[index, ]
          df$cell.pair <- factor(df$cell.pair, levels = unique(df$cell.pair))
          dotplot_fuc(df, high_freq_neglogpval_fishertest, out.dir, output_filename, genotype_colname)
        }
        } else if (dotplot_type == 3){
            ### Compare specific celltype vs all celltype
            if(is.null(celltype_a) | is.null(celltype_b)){
                print('Please assign celltype_a and celltype_b')
                stop()
            }
        
            df <- get(load(input_file_name))
            df1 <- df[grepl(celltype_a, df$cell.a, fixed=TRUE),]
            df1 <- df1[grepl(celltype_b, df1$cell.b, fixed=TRUE),]
            df <- df
            df$signif <- 0
            df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
            df$signif <- factor(df$signif, levels = c(0,1))
            # Select the minimum prop diff both side in mutant and wild type
            df = df[(df$mut_wt_prop_diff > mut_wt_prop_diff_max | df$mut_wt_prop_diff < mut_wt_prop_diff_min ),]

            # Select the iterecting pair which repeting minimum more than two time for each cell type pair
            #high_freq_int = names (table (df$interacting_pair)[table(df$interacting_pair) > 2])
            # Select the iterecting pair which repeting which are not repeting more than two time but have neglogpval_fishertest 1
            select_best_pvalue = df[df$mut_wt_neglogpval_fishertest > 0.99,]

            if (!is.null(top_ligand_receptor)) {
                # Select the top 10 rows based on 'mut_wt_neglogpval_fishertest'
                top_10_rows <- head(select_best_pvalue[order(select_best_pvalue$mut_wt_prop_diff, decreasing = TRUE), ], top_ligand_receptor)
                # Select the bottom 10 rows based on 'mut_wt_neglogpval_fishertest'
                bottom_10_rows <- tail(select_best_pvalue[order(select_best_pvalue$mut_wt_prop_diff, decreasing = TRUE), ], top_ligand_receptor)
                # Combine the top and bottom rows into a single DataFrame
                select_best_pvalue <- rbind(top_10_rows, bottom_10_rows)
                output_filename <- paste0("/dotplot_lig_rec_dt_",dotplot_type,"_top_",top_ligand_receptor,"_",celltype_a,"_",celltype_b,".pdf")
            } else {
              output_filename <- paste0("/dotplot_lig_rec_dt_",dotplot_type,"_all_sig_",celltype_a,"_",celltype_b,".pdf")
            }


            high_freq_neglogpval_fishertest = names (table (select_best_pvalue$interacting_pair))
            # select both list 
            df = df[df$interacting_pair %in% unique(c(high_freq_int,high_freq_neglogpval_fishertest)), ]
            dotplot_fuc(df, high_freq_neglogpval_fishertest, out.dir, output_filename, genotype_colname)

    } else if (dotplot_type == 4){
            ## Compare list of celltype with other cells
            df <- get(load(input_file_name))
            df <- df[grepl(paste(celltype_list,collapse="|"), df$cell.a) |
                grepl(paste(celltype_list,collapse="|"), df$cell.b), 
              ]
            df$signif <- 0
            df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
            df$signif <- factor(df$signif, levels = c(0,1))
            select_best_pvalue = df[df$mut_wt_neglogpval_fishertest > 0.99,]

            if (!is.null(top_ligand_receptor)) {
                # Select the top 10 rows based on 'mut_wt_neglogpval_fishertest'
                top_10_rows <- head(select_best_pvalue[order(select_best_pvalue$mut_wt_prop_diff, decreasing = TRUE), ], top_ligand_receptor)
                # Select the bottom 10 rows based on 'mut_wt_neglogpval_fishertest'
                bottom_10_rows <- tail(select_best_pvalue[order(select_best_pvalue$mut_wt_prop_diff, decreasing = TRUE), ], top_ligand_receptor)
                # Combine the top and bottom rows into a single DataFrame
                select_best_pvalue <- rbind(top_10_rows, bottom_10_rows)
                output_filename <- paste0("dotplot_lig_rec_dt_all_sig_only_top_",top_ligand_receptor,"_",Type_of_selcted_cells,"_vs_others.pdf")
            } else {
              output_filename <- paste0("dotplot_lig_rec_dt_all_sig_only_",Type_of_selcted_cells,"_vs_others.pdf")
            }

            high_freq_neglogpval_fishertest = names (table (select_best_pvalue$interacting_pair))
            df = df[df$interacting_pair %in% unique(c(high_freq_neglogpval_fishertest)), ]

            # Get the index that sorts the 'cell.a' column based on the presence of 'CAFs'
            index <- order(grepl(paste(celltype_list,collapse="|"), df$cell.a), decreasing = TRUE)
            # Reorder the DataFrame based on the index
            df <- df[index, ]
            df$cell.pair <- factor(df$cell.pair, levels = unique(df$cell.pair))
            dotplot_fuc(df, high_freq_neglogpval_fishertest, out.dir, output_filename, genotype_colname)
        } else if (dotplot_type == 5){
          for (gene_name in names(prog_list)){
              message(gene_name)
              genes_list = gene_list[[gene_name]]
              output_filename <- paste0("/dotplot_lig_rec_dt_",gene_name,"_all_sig_.pdf")
              df <- get(load(input_file_name))
              df$signif <- 0
              df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
              df$signif <- factor(df$signif, levels = c(0,1))
              select_best_pvalue = df[df$mut_wt_neglogpval_fishertest > 0.99,]
              
              ## Subset only selected genes
              subset_rows <- grepl(paste(genes_list, collapse = "|"), select_best_pvalue$interacting_pair)
              subset_df <- select_best_pvalue[subset_rows, ]
              select_best_pvalue <- subset_df
              
              high_freq_neglogpval_fishertest = names (table (select_best_pvalue$interacting_pair))
              df = df[df$interacting_pair %in% unique(c(high_freq_neglogpval_fishertest)), ]
              dotplot_fuc(df, high_freq_neglogpval_fishertest, out.dir, output_filename, genotype_colname)




              out.dir2 <- paste(out.dir, gene_name)
              dir.create(out.dir2)              
              df <- get(load(input_file_name))
              for(i in unlist(dimnames(table(df$cell.a)))){
                print(i)
                output_filename <- paste0("/dotplot_lig_rec_dt_",gene_name,"_all_sig_",i,".pdf")
                df <- get(load(input_file_name))
                df <- lig_rec_differential_highlvl_df_not_allzeros[
                  grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.a, fixed = TRUE) |
                  grepl(i, lig_rec_differential_highlvl_df_not_allzeros$cell.b, fixed = TRUE), 
                ]
                df$signif <- 0
                df$signif[df$mut_wt_neglogpval_fishertest > .99] <- 1
                df$signif <- factor(df$signif, levels = c(0,1))
                select_best_pvalue = df[df$mut_wt_neglogpval_fishertest > 0.99,]

                ## Subset only selected genes
                subset_rows <- grepl(paste(genes_list, collapse = "|"), select_best_pvalue$interacting_pair)
                subset_df <- select_best_pvalue[subset_rows, ]
                select_best_pvalue <- subset_df

                high_freq_neglogpval_fishertest = names (table (select_best_pvalue$interacting_pair))
                df = df[df$interacting_pair %in% unique(c(high_freq_neglogpval_fishertest)), ]
                
                # Get the index that sorts the 'cell.a' column based on the presence of specific celltype
                index <- order(grepl(i, df$cell.a, fixed = TRUE), decreasing = TRUE)
                # Reorder the DataFrame based on the index
                df <- df[index, ]
                df$cell.pair <- factor(df$cell.pair, levels = unique(df$cell.pair))
                dotplot_fuc(df, high_freq_neglogpval_fishertest, out.dir2, output_filename, genotype_colname)
              }

          }         
        }  else {
        print('Assign the dotplot_type between 1-4')
        print('Assign dotplot_type = 1, to generate all celltype combination ligand receptor pair dotplot')
        print('Assign dotplot_type = 2, to generate seperate ligand receptor pair dotplot for all celltype')
        print('Assign dotplot_type = 3, to generate ligand receptor pair dotplot between specific celltype')
        print('Assign dotplot_type = 4, to generate ligand receptor pair dotplot between one specific celltype vs others')
        print('Assign dotplot_type = 5, to generate ligand receptor pair dotplot from our genes list')
    }
}



dotplot_fuc <- function(df, high_freq_neglogpval_fishertest, out.dir, output_filename, genotype_colname){
    range <- max(abs(df$mut_wt_prop_diff))
    colors <- c(brewer_pal(palette = "Spectral", direction = -1)(7))
    pal <- gradient_n_pal(colors)
    custom_color_scale <- scale_fill_gradientn(
    colours = pal(c(0, rescale(seq_along(df$mut_wt_prop_diff)), 1)), # <- extra 0, 1 for out-of-bounds
    # limits = c(0, 6), breaks = 0:6,
    values = c(0,rescale(seq_along(df$mut_wt_prop_diff)),1), # <- extra 0, 1 again # also try changing to 7
    limits = c(-range,range)
    )

    if (length(c(high_freq_neglogpval_fishertest)) < 5){
      fig_length = 5
    } else if (length(c(high_freq_neglogpval_fishertest)) < 15){
      fig_length = 8
    } else{
      fig_length = length(c(high_freq_neglogpval_fishertest))/3
      if (fig_length < 8) {
          fig_length = 8
      }
    }
    

    if (length(table(df$cell.pair)) < 5){
        fig_width = 5
    } else if (length(table(df$cell.pair)) < 15){
        fig_width = 8
    } else {
        fig_width = length(table(df$cell.pair))/3
        if (fig_width < 8) {
            fig_width = 8
        }
    }

    

    colors <- c(brewer_pal(palette = "Spectral", direction = -1)(7))
    pal <- gradient_n_pal(colors)
    custom_color_scale <- scale_fill_gradientn(
    colours = pal(c(0, rescale(seq_along(df$mut_wt_prop_diff)), 1)), # <- extra 0, 1 for out-of-bounds
    # limits = c(0, 6), breaks = 0:6,
    values = c(0,rescale(seq_along(df$mut_wt_prop_diff)),1), # <- extra 0, 1 again # also try changing to 7
    limits = c(-range,range)
    )
    
    if (dim(df)[1] !=0){
    pdf (paste0(out.dir, output_filename), height = fig_length, width = fig_width)
    p <- ggplot(data = df, mapping = aes(x=cell.pair, y=interacting_pair, color=mut_wt_prop_diff, size=mut_wt_neglogpval_fishertest))+
    geom_point(shape = 21, aes(colour = as.factor(signif), fill = mut_wt_prop_diff))+
    scale_colour_manual(values=c("00FFFFFF", "black")) + theme_minimal() +
    custom_color_scale +
    theme(text = element_text(size=12), strip.text = element_text(size=22,face='bold'))+
    ggtitle(genotype_colname) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill = guide_colorbar(title = genotype_colname, labels = c("Red", "Blue"))) +
    labs(color = "Prop Diff", size = "Neg Log P-value")
    print(p)
    dev.off()
    }
}
