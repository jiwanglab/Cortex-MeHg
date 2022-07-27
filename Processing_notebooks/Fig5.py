import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
scv.settings.set_figure_params('scvelo')
#%%
adata = scv.read('Mix-SC.loom')
#%%
adata.var_names_make_unique()
### Control
#%%
metadata = pd.read_csv("Control_mix_metadata.csv")
genes = pd.read_csv("seurat_genes.csv")
cells = metadata['Barcode']
sample_adata = adata[cells.values, :]
sample_adata[:, np.intersect1d(genes.values,sample_adata.var_names.values)]
sample_adata.obsm['X_pca'] = np.column_stack([metadata['PC1'].values,
                                             metadata['PC2'].values])
sample_adata.obsm['X_umap'] = np.column_stack([metadata['UMAP1'].values,
                                               metadata['UMAP2'].values])

sample_adata.obs['seurat_clusters'] = metadata['seurat_clusters'].values

#%%

data = sample_adata.copy()
scv.pp.show_proportions(data)
scv.pp.filter_and_normalize(data, min_counts=20, min_counts_u=10, n_top_genes=3000)
scv.pp.moments(data)
scv.tl.velocity(data)
scv.tl.velocity_graph(data)
scv.tl.velocity_embedding(data, basis='umap')
scv.tl.velocity(data, mode='stochastic')


#%% Fig5a-b Control

scv.pl.velocity_embedding_stream(data, basis='umap', legend_loc='right margin',
                                     smooth=0.5, n_neighbors=20, density=1,
                                     color='seurat_clusters', size=30,groups=['Proliferative Radial Glial Cells','Radial Glial Cells','Intermediate Progenitors',
'Transitional Cells','Immature Neurons','Mature Neurons','Layer I Neurons','Interneurons','Endothelial Cells', 'Pericytes/Endothelial Cells'],
                                 palette={'Proliferative Radial Glial Cells':'#FD2D00','Radial Glial Cells':'#FD9900','Intermediate Progenitors':'#FDF900','Transitional Cells':'#04BC26', 'Immature Neurons':'#08E89D',
                                          'Mature Neurons':'#08E8DA', 'Layer I Neurons':'#08A8E8', 'Interneurons':'#D88BFF', 'Endothelial Cells':'#D700C0','Pericytes/Endothelial Cells':'#FC4BC9'})

### Velocity Fig5d Control
scv.tl.velocity_pseudotime(data)
scv.pl.scatter(data, color='velocity_pseudotime', cmap='gnuplot',save='Control_velocity.pdf')

### Fig5e Control
scv.pl.velocity(data, ["Pax6", "Tubb3","Tbr1"], color=['seurat_clusters'], save='Control_genes.pdf')

### Mercury

#%%

metadata = pd.read_csv("Mercury_mix_metadata.csv")
genes = pd.read_csv("seurat_genes.csv")
cells = metadata['Barcode']
sample_adata = adata[cells.values, :]
sample_adata[:, np.intersect1d(genes.values,sample_adata.var_names.values)]
sample_adata.obsm['X_pca'] = np.column_stack([metadata['PC1'].values,
                                             metadata['PC2'].values])
sample_adata.obsm['X_umap'] = np.column_stack([metadata['UMAP1'].values,
                                               metadata['UMAP2'].values])
sample_adata.obs['seurat_clusters'] = metadata['seurat_clusters'].values

adata.write_h5ad(filename='Mercury_Mix', compression=None, compression_opts=None, force_dense=None, as_dense=())

#%%

data = sample_adata
scv.pp.show_proportions(data)
scv.pp.filter_and_normalize(data, min_counts=20, min_counts_u=10, n_top_genes=3000)
scv.pp.moments(data)
scv.tl.velocity(data)
scv.tl.velocity_graph(data)
scv.tl.velocity_embedding(data, basis='umap')
scv.tl.velocity(data, mode='stochastic')

#%% Fig5a-b Mercury

scv.pl.velocity_embedding_stream(data, basis='umap', legend_loc='right margin',
                                     smooth=0.5, n_neighbors=20, density=1,
                                     color='seurat_clusters', size=30,
                                     dpi=300, groups=['Proliferative Radial Glial Cells','Radial Glial Cells','Intermediate Progenitors',
'Transitional Cells','Immature Neurons','Mature Neurons','Layer I Neurons','Interneurons','Endothelial Cells', 'Pericytes/Endothelial Cells'],
                                 palette={'Proliferative Radial Glial Cells':'#FD2D00','Radial Glial Cells':'#FD9900','Intermediate Progenitors':'#FDF900','Transitional Cells':'#04BC26', 'Immature Neurons':'#08E89D',
                                          'Mature Neurons':'#08E8DA', 'Layer I Neurons':'#08A8E8', 'Interneurons':'#D88BFF', 'Endothelial Cells':'#D700C0','Pericytes/Endothelial Cells':'#FC4BC9'},)


### Fig5d Mercury
scv.tl.velocity_pseudotime(data)
scv.pl.scatter(data, color='velocity_pseudotime', cmap='gnuplot',save='Mercury_velocity.pdf')

#%% Fig5f Mercury
scv.pl.velocity(data, ["Pax6", "Tubb3","Tbr1"], color=['seurat_clusters'], save='Mercury_genes.pdf')







import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
scv.settings.set_figure_params('scvelo')
#%%
adata = scv.read('Mix-SC.loom')
#%%
adata.var_names_make_unique()
### Control Cluster 1
#%%
metadata = pd.read_csv("Control_mix_metadata_cluster1.csv")
genes = pd.read_csv("seurat_genes_cluster1.csv")
cells = metadata['Barcode']
sample_adata = adata[cells.values, :]
sample_adata[:, np.intersect1d(genes.values,sample_adata.var_names.values)]
sample_adata.obsm['X_pca'] = np.column_stack([metadata['PC1'].values,
                                             metadata['PC2'].values])
sample_adata.obsm['X_umap'] = np.column_stack([metadata['UMAP1'].values,
                                               metadata['UMAP2'].values])

sample_adata.obs['seurat_clusters'] = metadata['seurat_clusters'].values

#%%

data = sample_adata.copy()
scv.pp.show_proportions(data)
scv.pp.filter_and_normalize(data, min_counts=10, min_counts_u=10, n_top_genes=3000)
scv.pp.moments(data)
scv.tl.velocity(data)
scv.tl.velocity_graph(data)
scv.tl.velocity_embedding(data, basis='umap')
scv.tl.velocity(data, mode='stochastic')

#%% Fig5g Control

scv.pl.velocity_embedding_stream(data, basis='umap', legend_loc='right margin',
                                     smooth=0.5, n_neighbors=20, density=1,
                                     color='seurat_clusters', size=30,
                                     dpi=300, groups=['Radial Glial Cells','Immature Neurons','Layer I Neurons','Pericytes'],
                                 palette={'Radial Glial Cells':'#F70D05','Immature Neurons':'#018916','Layer I Neurons':'#00D7E5','Pericytes':'#C133FA'})

### Fig5h Control
scv.tl.velocity_pseudotime(data)
scv.pl.scatter(data, color='velocity_pseudotime', cmap='gnuplot',save='Control1_velocity.pdf')

### Fig5i Control
scv.pl.velocity(data, ["Pax6","Tubb3", "Tbr1"],color=['seurat_clusters'],save='Control1_genes')

### Mercury Cluster 1

#%%

metadata = pd.read_csv("Mercury_mix_metadata_cluster1.csv")
genes = pd.read_csv("seurat_genes_cluster1.csv")
cells = metadata['Barcode']
sample_adata = adata[cells.values, :]
sample_adata[:, np.intersect1d(genes.values,sample_adata.var_names.values)]
sample_adata.obsm['X_pca'] = np.column_stack([metadata['PC1'].values,
                                             metadata['PC2'].values])
sample_adata.obsm['X_umap'] = np.column_stack([metadata['UMAP1'].values,
                                               metadata['UMAP2'].values])
sample_adata.obs['seurat_clusters'] = metadata['seurat_clusters'].values

#%%

data = sample_adata
scv.pp.show_proportions(data)
scv.pp.filter_and_normalize(data, min_counts=20, min_counts_u=10, n_top_genes=3000)
scv.pp.moments(data)
scv.tl.velocity(data)
scv.tl.velocity_graph(data)
scv.tl.velocity_embedding(data, basis='umap')
scv.tl.velocity(data, mode='stochastic')

### Fig5g Mercury
scv.pl.velocity_embedding_stream(data, basis='umap', legend_loc='right margin',
                                     smooth=0.5, n_neighbors=20, density=1,
                                     size=30, color='seurat_clusters',
                                     dpi=300, groups=['Radial Glial Cells','Immature Neurons','Layer I Neurons','Pericytes'],
                                 palette={'RGP1':'#F70D05','Immature Neurons':'#018916','Layer I Neurons':'#C133FA','Pericytes':'#00D7E5'}, save='mercury_cluster1.pdf')

### Fig5h Mercury
scv.tl.velocity_pseudotime(data)
scv.pl.scatter(data, color='velocity_pseudotime', cmap='gnuplot',save='Mercury1_velocity.pdf')

### Fig5j Mercury
scv.pl.velocity(data, ["Pax6","Tubb3","Tbr1"],color=['seurat_clusters'],save='Mercury1_genes')
#%%
