from getlines import distance_matrix
import sklearn.cluster

from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
import scipy.cluster.hierarchy as hcluster
import numpy as np

print(distance_matrix)

#clustering = SpectralClustering(n_clusters=3, assign_labels="discretize",random_state=0, affinity='precomputed').fit(distance_matrix)
#clustering1 = AgglomerativeClustering(n_clusters=None, compute_full_tree = True, distance_threshold = 100).fit(distance_matrix)
clustering2 = DBSCAN(eps=5, min_samples=2).fit(distance_matrix)


print(clustering2.labels_)






