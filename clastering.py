import sklearn.cluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs

import scipy.cluster.hierarchy as hcluster
import numpy as np
from math_funcs import CreateDistanceMatrix, HammingDistance
from getlines import writeinfile
from getlines import genome_array, genome_strings, distance_matrix

for genome in genome_array: # reworking the matrix
    if float(genome[2]) < 0.05:
        genome_array.remove(genome)

genome_strings = [genome[0] for genome in genome_array]
distance_matrix = CreateDistanceMatrix(genome_strings)
writeinfile(distance_matrix)

print(distance_matrix)

clustering = SpectralClustering(n_clusters=4, assign_labels="discretize",random_state=0, affinity='precomputed').fit(distance_matrix)
clustering1 = AgglomerativeClustering(n_clusters=None, compute_full_tree = True, distance_threshold = 10).fit(distance_matrix)
clustering2 = KMeans(n_clusters=5, random_state=1, n_init = 100).fit(distance_matrix)

labels = clustering2.labels_ # genome in genome_string[i] belongs to the claster in labels[i]
print(labels)

string_counter = 0
for string_counter in range(len(genome_array)): # here we append to genome_array, which cluster this genome belongs to
    genome_array[string_counter].append(labels[string_counter])
    string_counter = string_counter + 1

unique_labels = len(set(labels))
claster_sizes = [0] * unique_labels # we set a number of all labels - all but -1
for claster_number in labels:
    claster_sizes[claster_number] = claster_sizes[claster_number] + 1

X, y = make_blobs(n_samples=claster_sizes, centers=None, n_features=2, random_state=0, center_box = [-20, 20], cluster_std = 1) # creating a picture

plt.figure(figsize=(10, 10))

plt.scatter(X[:, 0], X[:, 1], c=y) # EACH POINT SHOULD CORRESPOND TO GENOME. No point corresponds to outliers

for point_counter in range(len(y)): #give to each non-(-1)-genome a corresponding point
    claster_num = y[point_counter]
    for genome in genome_array:
        if y[point_counter] == genome[3] and len(genome) == 4:
            genome.append(X[point_counter])
            break

for genome in genome_array:
    #plt.annotate(str(round(genome[2], 2)) + ' ' + str(genome[1]), (genome[4][0], genome[4][1])) # big annotate
    plt.annotate((round(genome[2], 2)), (genome[4][0], genome[4][1])) # small annotate

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














