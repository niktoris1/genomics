import sklearn.cluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
from sklearn.datasets import make_blobs
import scipy.cluster.hierarchy as hcluster
import numpy as np
from math_funcs import CreateDistanceMatrix
from getlines import writeinfile


from getlines import genome_array, genome_strings, distance_matrix



for genome in genome_array: # reworking the matrix
    if float(genome[2]) < 0.1:
        genome_array.remove(genome)

genome_strings = [genome[0] for genome in genome_array]
distance_matrix = CreateDistanceMatrix(genome_strings)
writeinfile(distance_matrix)

print(distance_matrix)

clustering = SpectralClustering(n_clusters=3, assign_labels="discretize",random_state=0, affinity='precomputed').fit(distance_matrix)
clustering1 = AgglomerativeClustering(n_clusters=None, compute_full_tree = True, distance_threshold = 100).fit(distance_matrix)
clustering2 = DBSCAN(eps=5, min_samples=2).fit(distance_matrix)

labels = clustering2.labels_ # genome in genome_string[i] belongs to the claster in labels[i]
print(labels)

unique_labels = len(set(labels))

if -1 in labels: #methods which assign some samples as outliers label them with -1. Here we get a true number of clasters
    unique_labels = unique_labels - 1

string_counter = 0
for genome in genome_array: # here we append to genome_array, which claster this genome belongs to
    genome.append(labels[string_counter])
    string_counter = string_counter + 1

claster_sizes = [0] * unique_labels
for claster_number in labels:
    claster_sizes[claster_number] = claster_sizes[claster_number] + 1


X, y = make_blobs(n_samples=claster_sizes, centers=None, n_features=2, random_state=0, center_box = [-20, 20], cluster_std = 2)

plt.scatter(X[:, 0], X[:, 1], c=y) # EACH POINT SHOLT CORRESPOND TO GENOME

for genome in genome_array:
    genome_cluster = genome[3]
    for counter in range(len(y)):
        if y[counter] == genome_cluster and len(genome) == 4:
            genome.append(X[counter]) # we also add to genome to coordinates of corresponding point X
            break

label_strings = [genome[0] for genome in genome_array]
label_shares = [genome[2] for genome in genome_array]

for i, genome in enumerate(genome_array):
    plt.annotate(round(label_shares[i], 2), (X[i][0], X[i][1]))

for first_genome in genome_array:
    for second_genome in genome_array:
        if first_genome is not second_genome and first_genome[1] == second_genome[1]:
            if first_genome[3] == second_genome[3]:
                color = 'green'
            else:
                color = 'red'
            x_values = [X[genome_array.index(first_genome), 0], X[genome_array.index(second_genome), 0]]
            y_values = [X[genome_array.index(first_genome), 1], X[genome_array.index(second_genome), 1]]
            plt.plot(x_values, y_values, color=color, linewidth=1)

figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
plt.show()














