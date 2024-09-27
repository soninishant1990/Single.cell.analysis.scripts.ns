
#Install Monocle3
# conda create --prefix /ahg/regevdata/projects/ICA_Lung/Nishant/conda/envs/monocle3 python=3.8
# conda install -c conda-forge r-base=4.1.2
# conda install -c bioconda r-monocle3
# install.packages('IRkernel')
# #install.packages('Seurat')

# install.packages(c("dplyr", "devtools", "ggpubr", "tidyverse", "RColorBrewer", "paletteer", "patchwork", "harmony"))
# conda install -c conda-forge r-seurat r-ggpubr r-gsa r-gdata r-tidyverse r-rcolorbrewer r-paletteer r-devtools r-nmf

# install.packages("remotes")
# remotes::install_github('satijalab/seurat-wrappers')
# install.packages('patchwork')

## https://ucdavis-bioinformatics-training.github.io/2020-August-Advanced-scRNAseq/data_analysis/adv_scrnaseq_monocle

## Run selected text in command line with enter+shift+a in visual studio
## "Ctrl+K, C" for selected line to comment out

#ssh -Y nsoni@login02.broadinstitute.org
ssh -Y nsoni@login01.broadinstitute.org
#screen -ls
#screen -r 'screen name' # if screen deteched
#screen -rd 'screen name' # if screen if atteched (or follow https://kb.iu.edu/d/ahrm)
#Ctrl + A and then Ctrl + D # to detached from screen
use UGER
ish -l os=RedHat7 -l h_vmem=16g -pe smp 2 -R y -binding linear:2
use UGER
conda activate monocle3
R
rm(list=ls()) # to remove Objects from an Environment



library (monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(ggplot2)
library(viridis)
library(dplyr)
library(paletteer)
library(harmony)

########################################################
########################################################
########################################################
# pediatric project path
## Microglia srt path
projDir = '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_3_pediatric/pediatric_histone_modification/_cellranger_filtered_Filter_400_1000_25/no_harmony/Myeloid_com_subset/sampleID_harmony/Myeloid_com_after_rev_doublet_subset/sampleID_harmony/Microglia_regress_out_subtype_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir, 'srt.rds'))

## MDM srt path
projDir = '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_3_pediatric/pediatric_histone_modification/_cellranger_filtered_Filter_400_1000_25/no_harmony/Myeloid_com_subset/sampleID_harmony/Myeloid_com_after_rev_doublet_subset/sampleID_harmony/MDM_com_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir, 'srt.rds'))


## Tumor srt path
projDir = '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_3_pediatric/pediatric_histone_modification/_cellranger_filtered_Filter_400_1000_25/no_harmony/Tumor_subset/sampleID_harmony/Tumor_after_rev_doublet_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir, 'srt.rds'))
srt$celltype = 'Tumor'
srt$sub_celltype = 'Tumor'


## Provide the column name for plot learn graph
metaGroupNames=c('Genotype_DM', 'sampleID', 'celltype', 'sub_celltype')

metaGroupNames=c('Genotype_DM', 'sampleID', 'celltype', 'sub_celltype')
reductionName = "sampleID_harmony_umap"


########################################################
########################################################
########################################################
## Genotype project
## Myeloid compartment
projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Myeloid.comp/Myeloid.comp_subset/sampleID_harmony/Myeloid.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir, 'srt.rds'))
projDirS = paste0 (projDir, 'Monocle1/')
dir.create (paste0(projDirS,'Plots/'), recursive=T, showWarnings = FALSE)

## Provide the column name for plot learn graph
metaGroupNames=c('groupID', 'sampleID', 'celltype', 'celltype')
reductionName = "sampleID_harmony_umap"

##########################################################
##########################################################
##########################################################
## MG 
projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Myeloid.comp/Myeloid.comp_subset/sampleID_harmony/Myeloid.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/MG.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir,'srt.rds'))
## Provide the column name for plot learn graph
metaGroupNames=c('groupID', 'sampleID', 'celltype', 'sub_celltype')
org = 'Mouse' ## 'Mouse' or 'Human'
reductionName = "sampleID_harmony_umap"
Monocle.dir.name = 'Monocle.with.harmony'

source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name)


### Run without harmony
srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
srt <- ScaleData(srt, features = VariableFeatures(object = srt))
srt <- RunPCA(srt, features = VariableFeatures(object = srt))
batch = 'no'

if (batch == 'no')
  {
  reductionName = 'umap'
  reductionSave = 'pca'
  reductionGraph = 'RNA_snn'
  } else {
  reductionSave = paste0(batch,'_harmony')
  reductionKey = paste0(batch,'harmonyUMAP_')
  reductionName = paste0 (batch,'_harmony_umap')
  reductionGraph = paste0 (batch,'_harmony_snn')
  }
sigPCs = 20
if (batch == 'no')
{
    srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)
} else {
    # Run Harmony
    srt = srt %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }
Monocle.dir.name = 'Monocle.without.harmony'

source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name,
org = 'Mouse')


##########################################################
##########################################################
##########################################################
### Disease associated MG
projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Myeloid.comp/Myeloid.comp_subset/sampleID_harmony/Myeloid.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/MG.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/Plots/Microglia.Disease.Associated/'
srt <- readRDS(paste0(projDir,'srt.rds'))
## Provide the column name for plot learn graph
metaGroupNames=c('groupID', 'sampleID', 'celltype', 'sub_celltype')
reductionName = "umap"

##########################################################
##########################################################
##########################################################
## MDM 
projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Myeloid.comp/Myeloid.comp_subset/sampleID_harmony/Myeloid.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/MDM.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir,'srt.rds'))
## Provide the column name for plot learn graph
metaGroupNames=c('groupID', 'sampleID', 'celltype', 'sub_celltype')
reductionName = "sampleID_harmony_umap"
srt$sub_celltype[srt$sampleID_harmony_snn_res.0.2 %in% c(0)] ='MDM.Disease.Associated'
srt$sub_celltype[srt$sampleID_harmony_snn_res.0.2 %in% c(1)] ='MDM.Pro.Inflammatory'
srt$sub_celltype[srt$sampleID_harmony_snn_res.0.2 %in% c(2)] ='MDM.Proliferating'

Monocle.dir.name = 'Monocle.with.harmony'
source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name)


### Run without harmony
srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
srt <- ScaleData(srt, features = VariableFeatures(object = srt))
srt <- RunPCA(srt, features = VariableFeatures(object = srt))
batch = 'no'

if (batch == 'no')
  {
  reductionName = 'umap'
  reductionSave = 'pca'
  reductionGraph = 'RNA_snn'
  } else {
  reductionSave = paste0(batch,'_harmony')
  reductionKey = paste0(batch,'harmonyUMAP_')
  reductionName = paste0 (batch,'_harmony_umap')
  reductionGraph = paste0 (batch,'_harmony_snn')
  }
sigPCs = 20
if (batch == 'no')
{
    srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)
} else {
    # Run Harmony
    srt = srt %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }
Monocle.dir.name = 'Monocle.without.harmony'

source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name)




##########################################################
##########################################################
##########################################################
## MDM and Monocytes
projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Myeloid.comp/Myeloid.comp_subset/sampleID_harmony/Myeloid.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir,'srt.rds'))
srt <- subset(srt, celltype == 'MDM' | celltype == 'Monocytes')
## Provide the column name for plot learn graph
metaGroupNames=c('groupID', 'sampleID', 'celltype', 'celltype')

### Run with harmony
srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
srt <- ScaleData(srt, features = VariableFeatures(object = srt))
srt <- RunPCA(srt, features = VariableFeatures(object = srt))
batch = 'sampleID'

if (batch == 'no')
  {
  reductionName = 'umap'
  reductionSave = 'pca'
  reductionGraph = 'RNA_snn'
  } else {
  reductionSave = paste0(batch,'_harmony')
  reductionKey = paste0(batch,'harmonyUMAP_')
  reductionName = paste0 (batch,'_harmony_umap')
  reductionGraph = paste0 (batch,'_harmony_snn')
  }
sigPCs = 20
if (batch == 'no')
{
    srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)
} else {
    # Run Harmony
    srt = srt %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }
Monocle.dir.name = 'MDM.Monocytes.Monocle.with.harmony'
source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name, org = 'Mouse')



### Run without harmony
srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
srt <- ScaleData(srt, features = VariableFeatures(object = srt))
srt <- RunPCA(srt, features = VariableFeatures(object = srt))
batch = 'no'

if (batch == 'no')
  {
  reductionName = 'umap'
  reductionSave = 'pca'
  reductionGraph = 'RNA_snn'
  } else {
  reductionSave = paste0(batch,'_harmony')
  reductionKey = paste0(batch,'harmonyUMAP_')
  reductionName = paste0 (batch,'_harmony_umap')
  reductionGraph = paste0 (batch,'_harmony_snn')
  }
sigPCs = 20
if (batch == 'no')
{
    srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)
} else {
    # Run Harmony
    srt = srt %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }
Monocle.dir.name = 'MDM.Monocytes.Monocle.without.harmony'
source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name, org = 'Mouse')




##########################################################
##########################################################
##########################################################
## MDM and Monocytes
projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Myeloid.comp/Myeloid.comp_subset/sampleID_harmony/Myeloid.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/Monocytes.MDM.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir,'srt.rds'))

## Provide the column name for plot learn graph
metaGroupNames=c('groupID', 'sampleID', 'celltype', 'sub_celltype')
Monocle.dir.name = 'Monocle.with.harmony'
projDirS = paste0 (projDir, Monocle.dir.name, '/')
dir.create (paste0(projDirS,'Plots/'), recursive=T, showWarnings = FALSE)


reductionName = "sampleID_harmony_umap"
source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name, org = 'Mouse')




### Run with harmony

srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
srt <- ScaleData(srt, features = VariableFeatures(object = srt))
srt <- RunPCA(srt, features = VariableFeatures(object = srt))
batch = 'sampleID'

if (batch == 'no')
  {
  reductionName = 'umap'
  reductionSave = 'pca'
  reductionGraph = 'RNA_snn'
  } else {
  reductionSave = paste0(batch,'_harmony')
  reductionKey = paste0(batch,'harmonyUMAP_')
  reductionName = paste0 (batch,'_harmony_umap')
  reductionGraph = paste0 (batch,'_harmony_snn')
  }
sigPCs = 15
if (batch == 'no')
{
    srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)
} else {
    # Run Harmony
    srt = srt %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }
Monocle.dir.name = 'Monocle.without.harmony'

source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name, org = 'Mouse')



##########################################################
##########################################################
##########################################################
## TUMOR
projDir <- '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/Tumor.comp/Tumor.comp_subset/sampleID_harmony/Tumor.Genotype.sample.NF1.EGFR.PDGFB_subset/sampleID_harmony/Tumor.after.removing.doublet_subset/sampleID_harmony/Tumor.after.regress.out.nfeature_subset/sampleID_harmony/'
srt <- readRDS(paste0(projDir,'srt.rds'))
## Provide the column name for plot learn graph
metaGroupNames=c('groupID', 'sampleID', 'celltype', 'celltype')
reductionName = "sampleID_harmony_umap"
Monocle.dir.name = 'Monocle.with.harmony'
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name)
# message ('Create dir for Monocle outputs')
# projDirS = paste0 (projDir, 'Monocle/')
# # dir.create (paste0(projDirS,'Plots/'), recursive=T, showWarnings = FALSE)

### Run without harmony

srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
srt <- ScaleData(srt, features = VariableFeatures(object = srt))
srt <- RunPCA(srt, features = VariableFeatures(object = srt))
batch = 'no'

if (batch == 'no')
  {
  reductionName = 'umap'
  reductionSave = 'pca'
  reductionGraph = 'RNA_snn'
  } else {
  reductionSave = paste0(batch,'_harmony')
  reductionKey = paste0(batch,'harmonyUMAP_')
  reductionName = paste0 (batch,'_harmony_umap')
  reductionGraph = paste0 (batch,'_harmony_snn')
  }
sigPCs = 20
if (batch == 'no')
{
    srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)
} else {
    # Run Harmony
    srt = srt %>% 
    RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
    RunUMAP (reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName, reduction.key=reductionKey)
    }
Monocle.dir.name = 'Monocle.without.harmony'

source('/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/LabCode/scrna_pipeline/Monocle3_func_ns.R')
Monocle.fun(srtS = srt, metaGroupNames = metaGroupNames, reductionName = reductionName, projDir = projDir, Monocle.dir.name = Monocle.dir.name, org = 'Mouse')











































### Monocle function
Monocle.fun <- function(srtS = NULL, metaGroupNames = NULL, reductionName = NULL, projDir = NULL, Monocle.dir.name = NULL){
    message ('Create dir for Monocle outputs')
    projDirS = paste0 (projDir, Monocle.dir.name, '/')
    dir.create (paste0(projDirS,'Plots/'), recursive=T, showWarnings = FALSE)

    # Check if any value from metaGroupNames is not present in colnames(srt@meta.data)
    if (any(!metaGroupNames %in% colnames(srtS@meta.data))) {
    stop("Code stopped: Some values from metaGroupNames are not present in colnames(srtS@meta.data).")
    }

    message ('convert seurat to cds obj')
    cds <- as.cell_data_set (srtS)

    ## Calculate size factors using built-in function in monocle3
    cds <- estimate_size_factors (cds)

    ## Add gene names into CDS
    cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srt[["RNA"]])
    cds <- cluster_cells(cds)

    ######################################################
    recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
    names(recreate.partitions) <- cds@colData@rownames
    recreate.partitions <- as.factor(recreate.partitions)
    #recreate.partitions

    cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
    #cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- srtS@reductions$umap@cell.embeddings
    cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- srtS@reductions[[reductionName]]@cell.embeddings

    message ('learn graph')
    cds <- learn_graph (cds, use_partition=FALSE, close_loop=FALSE)


    message ('Plot learn_graph')
    genes=NULL
    for(i in 1:length(metaGroupNames)){
        if(length(metaGroupNames[i]) != 1){
        pdf (paste0(projDirS, 'Plots/Monocle_srt_learn_graph_',dim (cds)[2],'_cells_umaps_',metaGroupNames[i],'.pdf'),8,5)
        print (plot_cells(cds,
                color_cells_by = metaGroupNames[i],
                label_principal_points = TRUE,
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=TRUE,
                graph_label_size=1.5,
                label_roots=TRUE,
                ))
        if (!is.null (genes)) print(plot_cells(cds,
                #color_cells_by = metaGroupName,
                label_cell_groups=FALSE,
                genes = genes,
                label_leaves=TRUE,
                label_branch_points=TRUE,
                graph_label_size=1.5))
        dev.off()
    }
    }


    message ('Plot show on show_trajectory_graph')
    p1 <- plot_cells(cds, show_trajectory_graph = TRUE, label_cell_groups=FALSE)
    pdf(paste0(projDirS , 'Plots/Monocle_srt_umaps_harmony.pdf'), width=8.5,height=7)
    print(wrap_plots(p1))
    dev.off()

    #####
    message ('Save cds object...')
    saveRDS (cds, paste0(projDirS, 'cds.rds'))
    cds = readRDS(paste0(projDirS, 'cds.rds'))


    #metaGroupNames=c('Genotype_DM', 'sampleID', 'celltype', 'sub_celltype')

    for(i in 1:length(metaGroupNames))
    {
        metaGroupName<-metaGroupNames[[i]]
        print(metaGroupName)
        pdf (paste0(projDirS,'Plots/UMAPs','_',metaGroupName,'_wHarmony_new_rigor.pdf'), width=8, height=5)
        umap=DimPlot(srt, reduction = reductionName, group.by = metaGroupName, label=FALSE,pt.size = .1,label.size = 5)
        plot(umap)
        dev.off()
    }



    ## Color cells by pseudotime
    message ('Plot show Color cells by pseudotime')
    cluster_levels <- as.numeric(cds@clusters$UMAP$clusters)
    max_cluster_level <- max(cluster_levels)
    root5 <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == max_cluster_level])) #4

    p1 <- plot_cells(root5,
            color_cells_by = "pseudotime",
            group_cells_by = "cluster",
            label_cell_groups = FALSE,
            label_groups_by_cluster=FALSE,
            label_leaves=FALSE,
            label_branch_points=FALSE,
            label_roots = FALSE,
            trajectory_graph_color = "black")

    pdf(paste0(projDirS , 'Plots/Color_cells_by_pseudotime1.pdf'), width=7,height=5)
    print(wrap_plots(p1))
    dev.off()

    png(paste0(projDirS , 'Plots/Color_cells_by_pseudotime1.png'), width = 1700, height = 1200,  res = 300)
    print(wrap_plots(p1))
    dev.off()


    ## Identify genes that change as a function of pseudotime. Monocle’s graph_test() function detects genes that vary over a trajectory.
    message ('Identify genes that change as a function of pseudotime')
    cds_graph_test_results <- graph_test(cds,
                                        neighbor_graph = "knn",
                                        cores = 2)
    head(cds_graph_test_results)
    deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

    p1 <- plot_cells(cds,
            genes = head(deg_ids),
            show_trajectory_graph = FALSE,
            label_cell_groups = FALSE,
            label_leaves = FALSE)

    png(paste0(projDirS , 'Plots/Identify_genes_that_change_as_a_function_of_pseudotime.png'), width=4000,height=2500, res=300)
    print(wrap_plots(p1))
    dev.off()


    # cds1 <- reduce_dimension(cds)
    # cds1 <- learn_graph (cds1, use_partition=FALSE, close_loop=FALSE)
    # p1 <- plot_cells(cds1, label_groups_by_cluster=FALSE,  color_cells_by = metaGroupNames[3])

    # pdf(paste0(projDirS , 'Plots/Reduce_dimensionality_visualize.pdf'), width=9,height=7)
    # print(wrap_plots(p1))
    # dev.off()
}


# print('NormalizeData')
# srt<- NormalizeData(object = srt)
# print('FindVariableFeatures')
# srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 3000)
# print('ScaleData')
# srt=ScaleData(srt)
# print('RunPCA')
# srt<- RunPCA(srt, features = VariableFeatures(object = srt))
# print('RunUMAP')
# srt <- RunUMAP(srt, reduction = "pca", dim=1:15)
# print('RunHarmony')
# srt <- RunHarmony(srt, c("sampleID"),assay.use="RNA")
# print('RunUMAP')
# srt <- RunUMAP(srt, reduction = "harmony", dim=1:15)
# srt

## Save seurat file
#saveRDS(srt, file = paste0(projDirS, 'srt_monocle3.rds'))
## Read seurat file
#srt <- readRDS(paste0(projDirS, 'srt_monocle3.rds'))
srtS = srt
message ('convert seurat to cds obj')
cds <- as.cell_data_set (srtS)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors (cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srt[["RNA"]])

cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds <- reduce_dimension(cds, umap.fast_sgd=FALSE, cores=1)
cds <- cluster_cells(cds)

######################################################
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
#recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
#cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- srtS@reductions$umap@cell.embeddings
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- srtS@reductions[[reductionName]]@cell.embeddings

message ('learn graph')
cds <- learn_graph (cds, use_partition=FALSE, close_loop=FALSE)


message ('Plot learn_graph')
genes=NULL
for(i in 1:length(metaGroupNames)){
    if(length(metaGroupNames[i]) != 1){
pdf (paste0(projDirS, 'Plots/Monocle_srt_learn_graph_',dim (cds)[2],'_cells_umaps_',metaGroupNames[i],'.pdf'),8,5)
print (plot_cells(cds,
           color_cells_by = metaGroupNames[i],
           label_principal_points = TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           label_roots=TRUE,
           ))
if (!is.null (genes)) print(plot_cells(cds,
           #color_cells_by = metaGroupName,
           label_cell_groups=FALSE,
           genes = genes,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5))
dev.off()
}
}


message ('Plot show on show_trajectory_graph')
p1 <- plot_cells(cds, show_trajectory_graph = TRUE, label_cell_groups=FALSE)
pdf(paste0(projDirS , 'Plots/Monocle_srt_umaps_harmony.pdf'), width=8.5,height=7)
print(wrap_plots(p1))
dev.off()

#####
message ('Save cds object...')
saveRDS (cds, paste0(projDirS, 'cds.rds'))
cds = readRDS(paste0(projDirS, 'cds.rds'))


#metaGroupNames=c('Genotype_DM', 'sampleID', 'celltype', 'sub_celltype')

for(i in 1:length(metaGroupNames))
{
    metaGroupName<-metaGroupNames[[i]]
    print(metaGroupName)
    pdf (paste0(projDirS,'Plots/UMAPs','_',metaGroupName,'_wHarmony_new_rigor.pdf'), width=8, height=5)
    umap=DimPlot(srt, reduction = reductionName, group.by = metaGroupName, label=FALSE,pt.size = .1,label.size = 5)
    plot(umap)
    dev.off()
}



## Color cells by pseudotime
message ('Plot show Color cells by pseudotime')
cluster_levels <- as.numeric(cds@clusters$UMAP$clusters)
max_cluster_level <- max(cluster_levels)
root5 <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == max_cluster_level])) #4

p1 <- plot_cells(root5,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "black")

pdf(paste0(projDirS , 'Plots/Color_cells_by_pseudotime1.pdf'), width=7,height=5)
print(wrap_plots(p1))
dev.off()

png(paste0(projDirS , 'Plots/Color_cells_by_pseudotime1.png'), width = 1700, height = 1200,  res = 300)
print(wrap_plots(p1))
dev.off()



message ('Plot show Color cells by pseudotime')
cluster_levels <- as.numeric(cds@clusters$UMAP$clusters)
max_cluster_level <- min(cluster_levels)
root5 <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == max_cluster_level])) #4


p1 <- plot_cells(root5,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "black")

# pdf(paste0(projDirS , 'Plots/Color_cells_by_pseudotime1.pdf'), width=7,height=5)
# print(wrap_plots(p1))
# dev.off()

png(paste0(projDirS , 'Plots/Color_cells_by_pseudotime2.png'), width = 1700, height = 1200,  res = 300)
print(wrap_plots(p1))
dev.off()



######################################
######################################
#### To reorder the trajectory 
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="Microglia.Homeostatic"){
  cell_ids <- which(colData(cds)[, "sub_celltype"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
root6 <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

p1 <- plot_cells(root6,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "black")

png(paste0(projDirS , 'Plots/Color_cells_by_pseudotime_order.by.Microglia.Homeostatic.png'), width = 1700, height = 1200,  res = 300)
print(wrap_plots(p1))
dev.off()


root6 <- extract_general_graph_ordering(cds, )





# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="Microglia.Proliferating"){
  cell_ids <- which(colData(cds)[, "sub_celltype"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
root6 <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

p1 <- plot_cells(root6,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "black")

png(paste0(projDirS , 'Plots/Color_cells_by_pseudotime_order.by.Microglia.Proliferating.png'), width = 1700, height = 1200,  res = 300)
print(wrap_plots(p1))
dev.off()




### Assign root for trajectory

# a helper function to identify the root principal points:
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)

  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

node_ids = get_correct_root_state(cds,
                                      cell_phenotype =
                                        'sub_celltype', "Microglia.Homeostatic")


## Identify genes that change as a function of pseudotime. Monocle’s graph_test() function detects genes that vary over a trajectory.
message ('Identify genes that change as a function of pseudotime')
cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "knn",
                                     cores = 2)
head(cds_graph_test_results)
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

p1 <- plot_cells(cds,
           genes = head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)

png(paste0(projDirS , 'Plots/Identify_genes_that_change_as_a_function_of_pseudotime.png'), width=4000,height=2500, res=300)
print(wrap_plots(p1))
dev.off()





cds1 <- reduce_dimension(cds)
cds1 <- learn_graph (cds1, use_partition=FALSE, close_loop=FALSE)
p1 <- plot_cells(cds1, label_groups_by_cluster=FALSE,  color_cells_by = metaGroupNames[3])

pdf(paste0(projDirS , 'Plots/Reduce_dimensionality_visualize.pdf'), width=9,height=7)
print(wrap_plots(p1))
dev.off()





