## Variable Set Up
projDir =  # path where this code is saved
fastqMetadata =  # metadata tab file from fastq2loom.py
runName = #just a title for file save names
metadata =  #metadata csv file from rnaVelocitySetUp.R
metaCol = ["sampleID", "celltype", 'condition', 'RNA_snn_res.0.8'] #clustering metadata to add, these are the columns I used for example

## Loom File Set Up
# Where the loom files saved, probably in the fast file folder
loomFilesdir = 
loomFiles = []

# This code assumes the use of a fastq metadata file used to construct the loom files
# If that does not exist, manually add loom file names to a loomFiles list and skip to SCVELO Set Up
import csv
sampleIDs = []

## ==================== ##
##  Using module 'csv'  ##
## ==================== ##
with open(fastqMetadata) as to_read:
    reader = csv.reader(to_read, delimiter = "\t")
    desired_column = [0] # SampleID column

    for row in reader:     # read one row at a time
        sampleIDs += list(row[i] for i in desired_column)   # build the output row (process)
        loomFiles.append(sampleIDs[-1] + ".loom")
                         
sampleIDs = sampleIDs[1:] # Skips header
sampleIDs

loomFiles = loomFiles[1:]
loomFiles # Should print something like ['KO1.loom', 'KO2.loom', 'KO3.loom', 'WT1.loom', 'WT2.loom', 'WT3.loom']

## SCVELO Set Up
import os
print (os.environ['CONDA_DEFAULT_ENV'])
print (os.getcwd())

import numpy
import pandas as pd
import scanpy as sc
numpy.version.version

import scvelo as scv
scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

## Load Data

## loompy used to convert fastqs to loom file for each sample
import loompy

if len(loomFiles) > 1:
    # Combine Samples: files variable set at top of page
    outputLoomFile = "combined.loom"
    os.chdir(loomFilesdir)
    if os.path.exists(outputLoomFile):
        os.remove(outputLoomFile) #helps overwrite the file to re-run
    loompy.combine(loomFiles, outputLoomFile, key="Accession")
else:
    # If there is only one loom file then just grab the name of the first thing in the list
    outputLoomFile = loomFiles[0]

# Load in loom file to Python/anndata (can also convert data from Seurat to Python / anndata)
# outputLoomFile = "combined.loom"
adata = scv.read(loomFilesdir + outputLoomFile, cache=True)
adata

# Good to check for error, if unspliced prop is 0, something went wrong and you have to recreate the loom files 
scv.pl.proportions(adata)


## Add in Metadata
# Load in metadata
# At this point, make sure you have run the rnaVelocitySetUp script 
os.chdir(projDir)
cell_meta = pd.read_csv(projDir + metadata, index_col='barcode')
cell_meta
cell_meta[metaCol]

# Adds umapData by comparing index column names in both
# UMAP embedding variables will always be the same
adata.obs['umap_1'] = cell_meta['umap_1']
adata.obs['umap_2'] = cell_meta['umap_2']
adata.obs

# Adds metadata by comparing index column names in both
# There will be NANs if the cell_metadata is only for a subset of data
adata.obs[metaCol] = cell_meta[metaCol].astype("string") #in case the clustering is read as numerical
adata.obs

#One way to remove NaN cells
adata = adata[adata.obs['umap_1'] == adata.obs['umap_1']] 
adata = adata[adata.obs['umap_2'] == adata.obs['umap_2']] 
adata.obs

# set umap embedding
adata.obsm['X_umap'] = numpy.vstack((adata.obs['umap_1'].to_numpy(), adata.obs['umap_2'].to_numpy())).T

# Load in metadata
## cell barcodes must match
os.chdir(projDir)
cell_meta = pd.read_csv(projDir + metadata, index_col='barcode')

# Adds umapData by comparing index column names in both
# UMAP embedding variables will always be the same
adata.obs['umap_1'] = cell_meta['umap_1']
adata.obs['umap_2'] = cell_meta['umap_2']

# Adds metadata by comparing index column names in both
adata.obs[metaCol] = cell_meta[metaCol].astype("string") #in case the clustering is read as numerical

# There will be NANs if the cell_metadata is only for a subset of data
# One way to remove NaN cells
adata = adata[adata.obs['umap_1'] == adata.obs['umap_1']] 
adata = adata[adata.obs['umap_2'] == adata.obs['umap_2']] 

# set umap embedding
adata.obsm['X_umap'] = numpy.vstack((adata.obs['umap_1'].to_numpy(), adata.obs['umap_2'].to_numpy())).T

# Plot a UMAP as sanity check:
sc.pl.umap(adata, color="sampleID", save = False)

# Plot a UMAP as sanity check:
for i in metaCol: # If there is more than one metagrouping
    sc.pl.umap(adata, color=[i], save = False) #save = "_" + runName + "_" + i + '_UMAP.png'


## RNA Velocity

# Estimate RNA Velocity
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# Project RNA Velocity 
# Projected onto any embedding on cellular level
scv.pl.velocity_embedding_stream(adata, basis='umap')
#pip install pandas==1.3.5

# Goal is to loop through and plot for each clustering added in the metadata
# Saves these figures in figures folder in the projDir, can set Save=False
os.chdir(projDir)
for i in metaCol: # If there is more than one metagrouping
    print(i)
    scv.pl.velocity_embedding_stream(adata, basis='umap', color=[i], dpi = 250, legend_loc = "right margin", save=runName + "_" + i +'_labelOff.png')


## Other Stats: Could be interesting, depending on what is being asked

## genes have cluster-specific differential velocity expression
# siginificantly higher/lower compared to the remaining population
# support the directionality and drive the velocities
scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

# Speed and coherence
# Speed = rate of differentiation (given by the length of the velocity vector)
# Coherence = confidence of the vector field
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save= runName + '_length_confidence.png')

# Can set a root here if asked of you 
# df = pd.DataFrame(adata.obs, index=adata.obs_names)
# root = adata.obs_names[adata.obs["RNA_snn_res.0.8"] == '0'][0]

# Velocity Graph and Velocity Pseudotime
# Average number of steps it takes to reach a cell after walking along the graph starting from the root cells
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save= runName + '_pseudotime.png')

# PAGA
# required to install igraph, if not done yet:
# !pip install python-igraph --upgrade --quiet

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups=metaCol[0])
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save= runName + '_PAGA.png')






