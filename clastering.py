import sklearn.cluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.cluster import MeanShift
from sklearn.cluster import AffinityPropagation



from sklearn.datasets import make_blobs

from External_data.get_external_data import blue_cluster, green_cluster, orange_cluster, red_cluster, vreden_cluster, white_cluster, vreden_small_cluster

import scipy.cluster.hierarchy as hcluster
import numpy as np
from math_funcs import CreateDistanceMatrix, HammingDistance
from getlines import writeinfile
from getlines import genome_array, genome_strings, distance_matrix

for genome in genome_array: # reworking the matrix
    if float(genome[2]) < 0.20:
        genome_array.remove(genome)

genome_strings = [genome[0] for genome in genome_array]
distance_matrix = CreateDistanceMatrix(genome_strings)
writeinfile(distance_matrix)

print(distance_matrix)

clustering = SpectralClustering(n_clusters=5, assign_labels="discretize",random_state=0, affinity='precomputed').fit(distance_matrix)
clustering1 = AgglomerativeClustering(n_clusters= None, compute_full_tree = True, affinity='precomputed', linkage = 'complete', distance_threshold=5).fit(distance_matrix)
clustering2 = KMeans(n_clusters=5, random_state=1, n_init = 100).fit(distance_matrix)
clustering3 = MeanShift().fit(distance_matrix)
clustering4 = AffinityPropagation().fit(distance_matrix)

labels = clustering1.labels_ # genome in genome_string[i] belongs to the cluster in labels[i]
print(labels)

string_counter = 0
for string_counter in range(len(genome_array)): # here we append to genome_array, which cluster this genome belongs to
    genome_array[string_counter].append(labels[string_counter])
    string_counter = string_counter + 1

unique_labels = len(set(labels))
cluster_sizes = [0] * unique_labels # we set a number of all labels - all but -1
for cluster_number in labels:
    cluster_sizes[cluster_number] = cluster_sizes[cluster_number] + 1

X, y = make_blobs(n_samples=cluster_sizes, centers=None, n_features=2, random_state=0, center_box = [-20, 20], cluster_std = 1) # creating a picture

plt.figure(figsize=(10, 10))

for point_counter in range(len(y)): #give to each non-(-1)-genome a corresponding point
    cluster_num = y[point_counter]
    for genome in genome_array:
        if y[point_counter] == genome[3] and len(genome) == 4:
            genome.append(X[point_counter])
            break

for genome in genome_array:
    if int(genome[1]) in blue_cluster:
        genome.append('blue')
    if int(genome[1]) in green_cluster:
        genome.append('green')
    if int(genome[1]) in orange_cluster:
        genome.append('orange')
    if int(genome[1]) in red_cluster:
        genome.append('red')
    if int(genome[1]) in white_cluster:
        genome.append('white')

for genome in genome_array:
    if int(genome[1]) in vreden_small_cluster:
        genome.append('vreden_small')

for genome in genome_array:
    plt.scatter(genome[4][0], genome[4][1], c=genome[5]) # EACH POINT SHOULD CORRESPOND TO GENOME. No point corresponds to outliers

for genome in genome_array:
    if len(genome) == 7:
        plt.annotate(str(round(genome[2], 2)) + 'V', (genome[4][0], genome[4][1])) # small annotate
    else:
        plt.annotate(str(round(genome[2], 2)), (genome[4][0], genome[4][1]))  # small annotate


for first_genome in genome_array:
    for second_genome in genome_array:
        if first_genome is not second_genome and first_genome[1] == second_genome[1]:
            if first_genome[3] == second_genome[3]:
                color = 'green'
            else:
                color = 'red'
            x_values = [first_genome[4][0], second_genome[4][0]]
            y_values = [first_genome[4][1], second_genome[4][1]]
            plt.annotate(HammingDistance(first_genome[0], second_genome[0]), ((first_genome[4][0] + second_genome[4][0]) / 2, (first_genome[4][1] + second_genome[4][1]) / 2))
            print(color + ' ' + str(HammingDistance(first_genome[0], second_genome[0])))
            plt.plot(x_values, y_values, color=color, linewidth=1)

plt.show()














