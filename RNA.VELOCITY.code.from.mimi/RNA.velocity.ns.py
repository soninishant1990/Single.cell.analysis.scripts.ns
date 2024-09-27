


#ssh -Y nsoni@login02.broadinstitute.org
ssh -Y nsoni@login01.broadinstitute.org
#screen -ls
#screen -r 'screen name' # if screen deteched
#screen -rd 'screen name' # if screen if atteched (or follow https://kb.iu.edu/d/ahrm)
#Ctrl + A and then Ctrl + D # to detached from screen
use UGER
ish -l os=RedHat7 -l h_vmem=64g -pe smp 2 -R y -binding linear:4
use UGER
conda activate rnavelocity
use Anaconda3
python
version('scvelo')


import csv
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


#### for Tumor Compaertment
## Variable Set Up
projDir =  '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/RNA.velocity/'# path where this code is saved
fastqMetadata = 'metadata.tab' # metadata tab file from fastq2loom.py
runName = 'Tumor.comprt.rna.velocity.csv' #just a title for file save names
metadata =  'Tumor.comprt.meta.file.csv' #metadata csv file from rnaVelocitySetUp.R
metaCol = ["sampleID", "celltype", 'groupID', 'sampleID_harmony_snn_res.0.4', 'Phase2'] #clustering metadata to add, these are the columns I used for example



#### for Microglia Compaertment
## Variable Set Up
projDir =  '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/RNA.velocity/'# path where this code is saved
fastqMetadata = 'metadata.tab' # metadata tab file from fastq2loom.py
runName = 'Microglia.comprt.rna.velocity.csv' #just a title for file save names
metadata =  'Microglia.comprt.meta.file.csv' #metadata csv file from rnaVelocitySetUp.R
metaCol = ["sampleID", "celltype", 'groupID', 'sub_celltype'] #clustering metadata to add, these are the columns I used for example


#### for MDM-Monocytes Compaertment
## Variable Set Up
projDir =  '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/RNA.velocity/'# path where this code is saved
fastqMetadata = 'metadata.tab' # metadata tab file from fastq2loom.py
runName = 'MDM.Monocytes.comprt.rna.velocity.csv' #just a title for file save names
metadata =  'MDM.Monocytes.comprt.meta.file.csv' #metadata csv file from rnaVelocitySetUp.R
metaCol = ["sampleID", "celltype", 'groupID', 'sub_celltype'] #clustering metadata to add, these are the columns I used for example



#### for Neutrophils Compaertment
## Variable Set Up
projDir =  '/ahg/regevdata/projects/ICA_Lung/Nishant/Dolores_prj_4_heterogeneity.between.genotype/RNA.velocity/'# path where this code is saved
fastqMetadata = 'metadata.tab' # metadata tab file from fastq2loom.py
runName = 'Neutrophils.comprt.rna.velocity.csv' #just a title for file save names
metadata =  'Neutrophils.comprt.meta.file.csv' #metadata csv file from rnaVelocitySetUp.R
metaCol = ["sampleID", 'groupID', 'sub_celltype'] #clustering metadata to add, these are the columns I used for example





## Loom File Set Up
# Where the loom files saved, probably in the fast file folder
loomFilesdir = projDir
file_names = os.listdir(loomFilesdir)
loomFiles = [file for file in file_names if file.endswith('.loom')]
print(loomFiles)
# len(loomFiles)
# # Get unique values using set
# unique_values = set(loomFiles)
# # Convert set back to list if needed
# unique_values_list = list(unique_values)
# print(unique_values_list)
# len(unique_values_list)

# Read the CSV file into a DataFrame
df = pd.read_csv(projDir+fastqMetadata, delimiter='\t')
print(df)
# Specify the desired column (SampleID column)
desired_column = 'name'

# Extract the values from the desired column
sampleIDs = df[desired_column].tolist()
sampleIDs

# len(sampleIDs)
# # Get unique values using set
# unique_values = set(sampleIDs)
# # Convert set back to list if needed
# unique_values_list = list(unique_values)
# print(unique_values_list)
# len(unique_values_list)


## loompy used to convert fastqs to loom file for each sample
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
outputLoomFile = "combined.loom"
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


# # Split each row name by underscore and extract the first value
# first_values = [row.split('_')[0] for row in missing_rownames]

# # Convert the list to a set to get unique values
# unique_first_values = set(first_values)

# print("Unique first values:", unique_first_values)

# import pandas as pd

# # Assuming df is your DataFrame
# search_value = '3920-19_AAACCCAGTCCAACGC'

# # Check if the search_value is present in any of the row names
# matches = adata.obs.index.str.contains(search_value)

# # Print the rows where the search_value is present in the row names
# print(adata.obs[matches].head())


# # Assuming df is your DataFrame
# adata.obs1 = adata.obs.drop_duplicates(inplace=True)
# print(adata.obs)

# print(adata.obs)
# print(cell_meta)

# # Assuming df is your DataFrame
# duplicates = adata.obs.index.duplicated(keep=False)

# # Check if there are any duplicate row names
# if duplicates.any():
#     print("Duplicate row names found.")
#     print("Duplicate row names:")
#     #print(adata.obs.index[duplicates])
# else:
#     print("No duplicate row names found.")

# Adds metadata by comparing index column names in both
# There will be NANs if the cell_metadata is only for a subset of data
adata.obs[metaCol] = cell_meta[metaCol].astype("string") #in case the clustering is read as numerical
adata.obs

#One way to remove NaN cells
adata = adata[adata.obs['umap_1'] == adata.obs['umap_1']] 
adata = adata[adata.obs['umap_2'] == adata.obs['umap_2']] 
adata.obs

# Remove duplicate row names from adata.obs DataFrame
#adata.obs.drop_duplicates(inplace=True)


# set umap embedding
adata.obsm['X_umap'] = numpy.vstack((adata.obs['umap_1'].to_numpy(), adata.obs['umap_2'].to_numpy())).T


# ### Preprocess the Data
# scv.pp.filter_genes(adata, min_shared_counts=20)
# scv.pp.normalize_per_cell(adata)
# scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
# scv.pp.log1p(adata)
# scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
# scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# scv.pp.remove_duplicate_cells(adata)

# # Load in metadata
# ## cell barcodes must match
# os.chdir(projDir)
# cell_meta = pd.read_csv(projDir + metadata, index_col='barcode')

# # Adds umapData by comparing index column names in both
# # UMAP embedding variables will always be the same
# adata.obs['umap_1'] = cell_meta['umap_1']
# adata.obs['umap_2'] = cell_meta['umap_2']

# # Adds metadata by comparing index column names in both
# adata.obs[metaCol] = cell_meta[metaCol].astype("string") #in case the clustering is read as numerical

# # There will be NANs if the cell_metadata is only for a subset of data
# # One way to remove NaN cells
# adata = adata[adata.obs['umap_1'] == adata.obs['umap_1']] 
# adata = adata[adata.obs['umap_2'] == adata.obs['umap_2']] 

# # set umap embedding
# adata.obsm['X_umap'] = numpy.vstack((adata.obs['umap_1'].to_numpy(), adata.obs['umap_2'].to_numpy())).T

# Plot a UMAP as sanity check:
#sc.pl.umap(adata, color="sampleID", save = False)
sc.pl.umap(adata, color="sampleID", save="umap_plot.png")

# Plot a# UMAP as sanity check:
for i in metaCol: # If there is more than one metagrouping
    sc.pl.umap(adata, color=[i], save = "_" + runName + "_" + i + '_UMAP.png') #save = "_" + runName + "_" + i + '_UMAP.png'


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
    scv.pl.velocity_embedding_stream(adata, basis='umap', color=[i], dpi = 300, legend_loc = "right margin", save=runName + "_" + i +'_labelOff.png')


### fOR TUMOR COMPARTMENT
cycling_color = "#006400"  # Dark green for cycling cells
non_cycling_color = "#ffc0cb"  # Light pink for non-cycling cells
cc_color = [cycling_color, non_cycling_color]
i = 'Phase2'
scv.pl.velocity_embedding_stream(adata, basis='umap', palette = cc_color, size=10, color=[i], dpi = 300, legend_loc = "right margin", save=runName + "_" + i +'_labelOff.png')

### For Microglia
#cycling_color = "#006400"  # Dark green for cycling cells
#non_cycling_color = "#ffc0cb"  # Light pink for non-cycling cells
cc_color = ['#AEC7E8FF', '#FFBB78FF', '#98DF8AFF', '#FF9896FF', '#C5B0D5FF']
i = 'sub_celltype'
scv.pl.velocity_embedding_stream(adata, basis='umap',  alpha=0.5, palette = cc_color, size=20, color=[i], dpi = 300, figsize =(7,7) , legend_loc = "right margin", save=runName + "_" + i +'_labelOff.png')


### For MDM/Monocytes
#cycling_color = "#006400"  # Dark green for cycling cells
#non_cycling_color = "#ffc0cb"  # Light pink for non-cycling cells
cc_color = ['#0055AAFF', '#C40003FF', '#00C19BFF', '#EAC862FF', '#7FD2FFFF', '#007ED3FF', '#B2DF8AFF', '#FFACAAFF']
i = 'sub_celltype'
scv.pl.velocity_embedding_stream(adata, basis='umap',  alpha=0.5, palette = cc_color, size=20, color=[i], dpi = 300, figsize =(7,7) , legend_loc = "right margin", save=runName + "_" + i +'_labelOff.png')


### For Neutrophils
cc_color = ['#BC3C29FF', '#0072B5FF', '#E18727FF', '#20854EFF', '#7876B1FF', '#6F99ADFF', '#FFDC91FF', '#EE4C97FF']
i = 'sub_celltype'
scv.pl.velocity_embedding_stream(adata, basis='umap',  alpha=0.5, palette = cc_color, size=20, color=[i], dpi = 300, figsize =(4,4) , legend_loc = "right margin", save=runName + "_" + i +'_labelOff.png')




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
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save= runName + '_length_confidence.png', dpi=300)

# Can set a root here if asked of you 
# df = pd.DataFrame(adata.obs, index=adata.obs_names)
# root = adata.obs_names[adata.obs["RNA_snn_res.0.8"] == '0'][0]

# Velocity Graph and Velocity Pseudotime
# Average number of steps it takes to reach a cell after walking along the graph starting from the root cells
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', size=5, cmap='gnuplot', save= runName + '_pseudotime.png', dpi=300)


# # Velocity graph and pseudotime
# scv.tl.velocity_pseudotime(adata)
# scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save= runName + '_velocity_pseudotime.png', dpi=300)


# scv.tl.score_genes_cell_cycle(adata)
# scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], , save= runName + '_cell_cycle.png')

# PAGA
# required to install igraph, if not done yet:
# !pip install python-igraph --upgrade --quiet

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups=metaCol[3])
scv.pl.paga(adata, basis='umap', size=5, alpha=.5, min_edge_width=2, node_size_scale=1.5, save= runName + '_PAGA.png', dpi=300)
scv.tl.paga(adata, groups=metaCol[4])
scv.pl.paga(adata, basis='umap', size=5, alpha=.5, min_edge_width=2, node_size_scale=1.5, save= runName + metaCol[4]+'_PAGA.png', dpi=300)

### For Neutrophils
scv.tl.paga(adata, groups=metaCol[2])
scv.pl.paga(adata, basis='umap', size=10, alpha=.8, min_edge_width=2, node_size_scale=1.5, save= runName + '_PAGA.png', dpi=300)


# Assume adata is your AnnData object
#scv.write(adata, "Tumur.comp.output.after.rna.velocity.cal.loom")
adata.write("Tumur.comp.output.after.rna.velocity.cal.loom")
#adata.write_loom("Tumur.comp.output.after.rna.velocity.cal.loom")
adata1 = scv.read(loomFilesdir + "Tumur.comp.output.after.rna.velocity.cal.loom", cache=True)
#adata = anndata.read_loom("Tumur.comp.output.after.rna.velocity.cal.loom")