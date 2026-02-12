# This code is based off of the code done in Vierstra et al 2020 (https://github.com/vierstralab/motif-clustering/blob/master/Workflow_v2.1beta-human.ipynb) 
# and uses the vertebrate HOMER motif file (turned into MEME format using the universalmotif package in R)

import pandas as pd
import numpy as np

# Load up the TOMTOM results (and remove the footer that isn't needed)
tomtom = pd.read_table('tomtom.homer.all.txt', skipfooter=3)
sim = tomtom.pivot_table(index='Query_ID', columns='Target_ID', values='E-value', fill_value=np.nan)

# Replace the occurances of inf (such as CTCF to itself) to the minimum value otherwise. Yeah
sim = sim.replace(0, sim[sim != 0].min(axis=None))

x = sim.values

w = np.triu(x) + np.triu(x, 1).T
v = np.tril(x) + np.tril(x, -1).T

sim.iloc[:,:] = np.nanmin(np.dstack([w, v]), axis=2)

### NOTE: The above returns a warning as it does in the Vierstra code
# <stdin>:1: RuntimeWarning: All-NaN slice encountered

sim.fillna(100, inplace=True)
sim = -np.log10(sim)
sim[np.isinf(sim)] = 10


# Cluster the square matrix
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram

Z = linkage(sim, method = 'complete', metric = 'correlation')

cl = fcluster(Z, 0.7, criterion='distance')
o = dendrogram(Z, no_plot=True)['leaves']

motif_annot_df = pd.DataFrame({'motif_id':sim.index, 'cluster':cl}).set_index('motif_id')

# With the clusters looking pretty decent, write the object to file! Woot woot
motif_annot_df.to_csv('homer_mot_cluster_python.txt', header=True, index=True, sep='\t')