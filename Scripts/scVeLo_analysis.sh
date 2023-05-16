#https://smorabit.github.io/tutorials/8_velocyto/
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("/home/c_x_he/CdCS/scVeLo/counts.mtx")

# create anndata object
adata = anndata.AnnData(X=X.transpose().tocsr())

# load cell metadata:
cell_meta = pd.read_csv("/home/c_x_he/CdCS/scVeLo/metadata.csv")

# load gene names:
with open("/home/c_x_he/CdCS/scVeLo/gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("/home/c_x_he/CdCS/scVeLo/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['orig.ident'], frameon=False, save=True)

# save dataset as anndata format
adata.write('/home/c_x_he/CdCS/scVeLo/seurat_to_anndata.h5ad')



################Step 0: Constructing spliced and unspliced counts matrices
velocyto run10x -m /home/c_x_he/CdCS/velocyto/repeat_msk.gtf /home/c_x_he/CdCS/cellranger_outs/CdCS/WES1377_1_d90 /home/Camila/opt/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf

#Generating output file /home/c_x_he/CdCS/cellranger_outs/CdCS/WES1377_1_d90/velocyto/WES1377_1_d90.loom

velocyto run10x -m /home/c_x_he/CdCS/velocyto/repeat_msk.gtf /home/c_x_he/CdCS/cellranger_outs/CdCS/WES1377_2_d90 /home/Camila/opt/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf


velocyto run10x -m /home/c_x_he/CdCS/velocyto/repeat_msk.gtf /home/c_x_he/CdCS/cellranger_outs/CdCS/WES1378_d90 /home/Camila/opt/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf

velocyto run10x -m /home/c_x_he/CdCS/velocyto/repeat_msk.gtf /home/c_x_he/CdCS/cellranger_outs/CdCS/K11824_d90 /home/Camila/opt/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf

###########Step 1: Load data
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

# reload dataset
adata = sc.read_h5ad('/home/c_x_he/CdCS/scVeLo/seurat_to_anndata.h5ad')

# load loom files for spliced/unspliced matrices for each sample:
# make variable names unique
ldata1 = scv.read("/home/c_x_he/CdCS/cellranger_outs/CdCS/WES1377_1_d90/velocyto/WES1377_1_d90.loom", cache=True)
ldata1.var_names_make_unique()

ldata2 = scv.read("/home/c_x_he/CdCS/cellranger_outs/CdCS/WES1377_2_d90/velocyto/WES1377_2_d90.loom", cache=True)
ldata2.var_names_make_unique()

ldata3 = scv.read("/home/c_x_he/CdCS/cellranger_outs/CdCS/K11824_d90/velocyto/K11824_d90.loom", cache=True)
ldata3.var_names_make_unique()

# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_9' for bc in barcodes]
ldata3.obs.index = barcodes

# concatenate the three loom
ldata = ldata1.concatenate([ldata2, ldata3])

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# plot umap to check
sc.pl.umap(adata, color='Cell.type', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')

####Part 2: Computing RNA velocity using scVelo
scv.pl.proportions(adata, groupby='Cell.type', save='_celltypes-table.pdf')

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)


#####Part 2.1: Visualize velocity fields
scv.pl.velocity_embedding_grid(adata, basis='umap', color='Cell.type', save='embedding_grid.pdf', title='CdCS RNA Velocity fields', scale=0.25)

scv.pl.velocity_embedding_stream(adata, basis='umap', color=['Cell.type', 'orig.ident'], save='embedding_stream.pdf', title='CdCS RNA Velocity fields')

#####Error when trying to save graph
scv.pl.velocity(adata, var_names=['SOX2'], color='Cell.type')

####Part 3: Downstream analysis
scv.tl.rank_velocity_genes(adata, groupby='Cell.type', min_corr=.3)

dataframe = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
dataframe.head()

#write.csv(dataframe, save="/home/c_x_he/CdCS/gene-velocity-rank.csv")

#view predicted projections by 2 celltypes 
kwargs = dict(frameon=False, size=10, linewidth=1.5,
              add_outline='AF6, AF1')


scv.pl.scatter(adata, dataframe['NPC_Ast+'][:3], ylabel='Ast', frameon=False, color='Cell.type', size=10, linewidth=1.5, save='ast-graph.pdf')

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save='Velocity-confidence.pdf')

#pseudotime
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save='Velocity-psudotime.pdf')


######Part 4: Analyzing a specific cell population

cur_celltypes = ['AF1', 'AF2', 'AF3', 'AF4', 'AF5', 'AF6', 'PF1', 'PF2']
adata_subset = adata[adata.obs['celltype'].isin(cur_celltypes)]
sc.pl.umap(adata_subset, color=['celltype', 'condition'], frameon=False, title=['', ''])

sc.pp.neighbors(adata_subset, n_neighbors=15, use_rep='X_pca')

# pre-process
scv.pp.filter_and_normalize(adata_subset)
scv.pp.moments(adata_subset)
