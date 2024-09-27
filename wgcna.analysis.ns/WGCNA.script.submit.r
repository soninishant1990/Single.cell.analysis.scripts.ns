# load libraries
message ('Load R packages')
packages = c('ggplot2','Seurat','dplyr','Matrix','ggpubr',
  'biomaRt','gplots','patchwork','ComplexHeatmap',
  'RColorBrewer','ggrepel','fgsea',
  'tidyr','gdata','GSA',
  'DropletUtils', 
  'harmony','scater',
  'destiny','plotly',
  'GO.db', 'org.Mm.eg.db','org.Hs.eg.db',
  'viridis',
  'BiocParallel',
  'ggridges',
  'clusterProfiler',
  'NMF',
  'WGCNA',
  'scWGCNA','igraph','tidyverse', 'paletteer'
)
lapply(packages, require, character.only = TRUE)

projDir = '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_3_pediatric/pediatric_histone_modification/_cellranger_filtered_Filter_400_1000_25/no_harmony/Myeloid_com_subset/sampleID_harmony/Myeloid_com_after_rev_doublet_subset/sampleID_harmony/'

srt <- readRDS(paste0(projDir, 'srt.rds'))

### Change path and the name of the dataset
## Create new path
out_dic <- projDir
dir <- paste0(out_dic,'wgcna_analysis_ns/')
dir.create(dir)

## srt file name
srt_obj_filename = paste0(projDir,"srt.rds")
## data organism # 'mouse' or 'human'
org='mouse' # organism
## reduction name and save into pca
reductionSave = 'pca'
reductionName = 'sampleID_harmony_umap'

## metgroup for WGCNA analysis - sampleid and groupid column
metacells_groups = c('Genotype_DM','sampleID_harmony_snn_res.0.2') # set metagroup in which to find metacells (usually your clustering / celltypes)
metaGroupNames = c('orig.ident','sub_celltype','Genotype_DM') # set of metagroups for generating the boxplots

### single cell analysis usefull script and WGCNA script load
useful_fun_path <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_3_pediatric/LabCode/useful_functions.R'
wgcna_source_code_path <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_3_pediatric/LabCode/scrna_pipeline/wgcna_analysis_ns.R'

## Assign parameter for job submit
mail_id <- 'nsoni@broadinstitute.org'
vmem <- '32g'
smp<- '2'
h_rt <- '24:00:00'
conda_env <- 'scrna'

## Generate R script for WGCNA SCRIPT
f <- file(paste0(dir,"wgcna_script.R"), open = "wt")
cat <- function(...){
    base::cat(..., file=f)
    }
cat(paste0('

library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(clusterProfiler)

myargs = commandArgs(trailingOnly=TRUE)
softPower=strtoi(myargs[1], base = 0L)
deepSplit_val=strtoi(myargs[2], base = 0L)
mergeCutHeight = strtoi(myargs[3], base = 0L)
metacells_k_val=strtoi(myargs[4], base = 0L)
max_shared_val=strtoi(myargs[5], base = 0L)
minModuleSize_val=strtoi(myargs[6], base = 0L)
nfeatures = strtoi(myargs[7], base = 0L)


wgcna_dir <- "',dir,'"
#dir.create(wgcna_dir)

srt <- readRDS(paste0("',srt_obj_filename,'"))
#projDir = projDir # directory where to save outputs 
srt_wgcna = srt # seurat object
org="',org,'" # organism
reductionName = "',reductionName,'" # non-linear dimensionality reduction to use
reductionSave = "',reductionSave,'" # linear dimensionality reduction to use
force=TRUE # force re-running WGCNA 
do.fgsea=TRUE # Run pathway enrichments
powerTable=FALSE # plot the powerTable plot. Can take a while to generate
do.plots=TRUE # generate plots
softPower=softPower # Set the softPower
deepSplit=deepSplit_val # Set this to increase decrease the number of modules idendified. 1-4
mergeCutHeight = 0.20 # Height below which two modules are merged
metacells_k = metacells_k_val # number of cells to create pseudobulks
max_shared = max_shared_val
minModuleSize=minModuleSize_val
metacells_groups = c("',metacells_groups[[1]],'","',metacells_groups[[2]],'") # set metagroup in which to find metacells (usually your clustering / celltypes)
metaGroupNames = c("',metaGroupNames[[1]],'","',metaGroupNames[[2]],'","',metaGroupNames[[3]],'") # set of metagroups for generating the boxplots
module_pal = viridis::viridis (length(unique (srt_wgcna$orig.ident)))
genes.keep = VariableFeatures (FindVariableFeatures (srt_wgcna, nfeatures = nfeatures)) # Genes to use to compute WGCNA
#genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeatures=length(rownames(srt))))
enricher_universe = rownames (srt_wgcna) # Genes to use as background for pathway enrichments analysis
message("parameter loaded")
source ("',useful_fun_path,'")
source("',wgcna_source_code_path,'")'))
close(f)

## Generate sh file for job submitting
 f <- file(paste0(dir,"/wgcna_job.sh"), open = "wt")
cat <- function(...){
    base::cat(..., file=f)
    }
cat(paste0('#!/bin/bash
#$ -M ',mail_id,'
#$ -N job.wgcna
#$ -cwd
#$ -q broad
#$ -l h_vmem=',vmem,'
#$ -pe smp ',smp,'
#$ -binding linear:2
#$ -l h_rt=',h_rt,'
#$ -e cpdb.err
#$ -o ',dir,'

source /broad/software/scripts/useuse
use .anaconda3-2022.10
source activate ',conda_env,'

softPower=${1}
deepSplit_val=${2}
mergeCutHeight = ${3}
metacells_k_val=${4}
max_shared_val=${5}
minModuleSize_val=${6}
nfeatures = ${7}

Rscript ',dir,'wgcna_analysis_ns.R ${1} ${2} ${3} ${4} ${5} ${6} ${7}
'))
close(f)

## Change the permission for sh file
system(paste0('chmod u+x ',dir,'wgcna_job.sh'), wait=FALSE)



## Submit job with multiple parameter
#####
softPower=c(12)
deepSplit_val=c(4)
mergeCutHeight = c(0.20)
metacells_k_val=c(40)
max_shared_val=c(25)
minModuleSize_val=c(50)
nfeatures <- c(10000)

# Loop to process the first values
for (i in 1:length(softPower)) {
  current_softPower <- softPower[i]
  current_deepSplit_val <- deepSplit_val[i]
  current_mergeCutHeight <- mergeCutHeight[i]
  current_metacells_k_val <- metacells_k_val[i]
  current_max_shared_val <- max_shared_val[i]
  current_minModuleSize_val <- minModuleSize_val[i]
  current_nfeatures <- nfeatures[i]
    system(paste0('qsub ',dir,'wgcna_job.sh ',current_softPower,' ',current_deepSplit_val,' ',current_mergeCutHeight,' ',
             current_metacells_k_val,' ',current_max_shared_val,' ',current_minModuleSize_val,' ',current_nfeatures), wait=FALSE)
}

#####
softPower=c(12, 12, 12, 12, 12, 12, 12)
deepSplit_val=c(4, 4, 4, 4, 4, 4, 4)
mergeCutHeight = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20)
metacells_k_val=c(40, 40, 40, 40, 40, 40, 40)
max_shared_val=c(25, 25, 25, 25, 25, 25, 25)
minModuleSize_val=c(15, 20, 25, 30, 40, 45, 50)
nfeatures <- c(10000, 10000, 10000, 10000, 10000, 10000, 10000)

# Loop to process the first values
for (i in 1:length(softPower)) {
  current_softPower <- softPower[i]
  current_deepSplit_val <- deepSplit_val[i]
  current_mergeCutHeight <- mergeCutHeight[i]
  current_metacells_k_val <- metacells_k_val[i]
  current_max_shared_val <- max_shared_val[i]
  current_minModuleSize_val <- minModuleSize_val[i]
  current_nfeatures <- nfeatures[i]
    system(paste0('qsub ',dir,'wgcna_job.sh ',current_softPower,' ',current_deepSplit_val,' ',current_mergeCutHeight,' ',
             current_metacells_k_val,' ',current_max_shared_val,' ',current_minModuleSize_val,' ',current_nfeatures), wait=FALSE)
}
#####
softPower=c(12, 12, 12, 12, 12, 12, 12)
deepSplit_val=c(4, 4, 4, 4, 4, 4, 4)
mergeCutHeight = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20)
metacells_k_val=c(50, 50, 50, 50, 50, 50, 50)
max_shared_val=c(25, 25, 25, 25, 25, 25, 25)
minModuleSize_val=c(15, 20, 25, 30, 40, 45, 50)
nfeatures <- c(10000, 10000, 10000, 10000, 10000, 10000, 10000)

# Loop to process the first values
for (i in 1:length(softPower)) {
  current_softPower <- softPower[i]
  current_deepSplit_val <- deepSplit_val[i]
  current_mergeCutHeight <- mergeCutHeight[i]
  current_metacells_k_val <- metacells_k_val[i]
  current_max_shared_val <- max_shared_val[i]
  current_minModuleSize_val <- minModuleSize_val[i]
  current_nfeatures <- nfeatures[i]
    system(paste0('qsub ',dir,'wgcna_job.sh ',current_softPower,' ',current_deepSplit_val,' ',current_mergeCutHeight,' ',
             current_metacells_k_val,' ',current_max_shared_val,' ',current_minModuleSize_val,' ',current_nfeatures), wait=FALSE)
}
#####
softPower=c(12, 12, 12, 12, 12, 12, 12)
deepSplit_val=c(4, 4, 4, 4, 4, 4, 4)
mergeCutHeight = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20)
metacells_k_val=c(60, 60, 60, 60, 60, 60, 60)
max_shared_val=c(25, 25, 25, 25, 25, 25, 25)
minModuleSize_val=c(15, 20, 25, 30, 40, 45, 50)
nfeatures <- c(10000, 10000, 10000, 10000, 10000, 10000, 10000)

# Loop to process the first values
for (i in 1:length(softPower)) {
  current_softPower <- softPower[i]
  current_deepSplit_val <- deepSplit_val[i]
  current_mergeCutHeight <- mergeCutHeight[i]
  current_metacells_k_val <- metacells_k_val[i]
  current_max_shared_val <- max_shared_val[i]
  current_minModuleSize_val <- minModuleSize_val[i]
  current_nfeatures <- nfeatures[i]
    system(paste0('qsub ',dir,'wgcna_job.sh ',current_softPower,' ',current_deepSplit_val,' ',current_mergeCutHeight,' ',
             current_metacells_k_val,' ',current_max_shared_val,' ',current_minModuleSize_val,' ',current_nfeatures), wait=FALSE)
}






### change module name Regenerate boxplot

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

projDirW =  '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_3_pediatric/pediatric_histone_modification/_cellranger_filtered_Filter_400_1000_25/no_harmony/Myeloid_com_subset/sampleID_harmony/Myeloid_com_after_rev_doublet_subset/sampleID_harmony/MDM_com_subset/sampleID_harmony/wgcna_analysis_ns/WGCNA_sp_12_ds_4_mch_0.2max_shared_25_minModuleSize_50_metasize_60_genes_10000/'

srt_wgcna <- readRDS(paste0(projDirW, 'srt_wgcna.rds'))
modules <- GetModules(srt_wgcna)
modulesL = split (modules$gene_name, modules$module)

new_colnames = c('EMT', 
                'gray',
                'Hypoxia_Glycolysis',
                'G2.M',
                'IL2_STAT5', 
                'TNFA_Signaling',
                'Unknown1',
                'Interferon',
                'Unknown2', 
                'G1.S', 
                'MTROC1',
                'Unknown3', 
                'MYC_target',
                'Neurogenesis',
                'KRAS_Signaling')

colnames(srt_wgcna@meta.data)[46:62] <- new_colnames

metaGroupNames = c('orig.ident','Genotype_DM')
meta_modules_names = new_colnames # remove grey module (unassigned)
umap_df = data.frame (srt_wgcna[[reductionName]]@cell.embeddings, srt_wgcna@meta.data[,meta_modules_names])
umap_p1 = lapply (meta_modules_names, function(x) ggplot(data = umap_df) + 
geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
theme_bw() + ggtitle (x))

### Generate boxplots per meta groups
ccomp_df = srt_wgcna@meta.data[,new_colnames, drop=FALSE]
#ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 
ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar
#color_bar <- paletteer_c("grDevices::rainbow", dim(table(srt_wgcna@meta.data[,metaGroupNames[2]])))

color_bar <- c("BS.Histone.H3.1K27M" = "purple",
                   "BS.Histone.H3.3K27M" = "blue",
                   "BS.Histone.H3.wt" = "red",
                   "Histone.H3.wt" = "grey")
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

png (paste0 (projDirW,'Plots/WGCNA_module_scores_umap_new_ann_at_genotype.png'), width = 9000, height = 4000, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))) / 
wrap_plots (box_p, ncol= ifelse (length(box_p) > 8,ceiling(length(box_p)/2),length(box_p))) + plot_layout(heights=c(1,1.5)))
dev.off()

png (paste0 (projDirW,'Plots/WGCNA_module_scores_boxplot_new_ann_at_genotype.png'), width = 5700, height = 2800, pointsize=10, res = 300, type="cairo")
print (wrap_plots (box_p, ncol= ifelse (length(box_p) > 8,ceiling(length(box_p)/2),length(box_p))))
dev.off()

metaGroupNames = c('orig.ident','sub_celltype')
meta_modules_names = new_colnames # remove grey module (unassigned)
umap_df = data.frame (srt_wgcna[[reductionName]]@cell.embeddings, srt_wgcna@meta.data[,meta_modules_names])
umap_p1 = lapply (meta_modules_names, function(x) ggplot(data = umap_df) + 
geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
theme_bw() + ggtitle (x))

### Generate boxplots per meta groups
ccomp_df = srt_wgcna@meta.data[,new_colnames, drop=FALSE]
#ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 
ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)

## Get the color bar
color_bar <- paletteer_c("grDevices::rainbow", dim(table(srt_wgcna@meta.data[,metaGroupNames[2]])))

#color_bar <- c("BS.Histone.H3.1K27M" = "purple",
#                   "BS.Histone.H3.3K27M" = "blue",
#                   "BS.Histone.H3.wt" = "red",
#                   "Histone.H3.wt" = "grey")


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

png (paste0 (projDirW,'Plots/WGCNA_module_scores_umap_new_ann_at_sub_celltype.png'), width = 9000, height = 4000, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))) / 
wrap_plots (box_p, ncol= ifelse (length(box_p) > 8,ceiling(length(box_p)/2),length(box_p))) + plot_layout(heights=c(1,1.5)))
dev.off()

png (paste0 (projDirW,'Plots/WGCNA_module_scores_boxplot_new_ann_at_sub_celltype.png'), width = 5700, height = 3000, pointsize=10, res = 300, type="cairo")
print (wrap_plots (box_p, ncol= ifelse (length(box_p) > 8,ceiling(length(box_p)/2),length(box_p))))
dev.off()























