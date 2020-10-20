from math_funcs import CreateDistanceMatrix

def getgenomes():

    file = open("resulting_genomes.txt", "r")

    raw_genomes = [line.split() for line in file]

    genomes = [] #format is [genome, sample, share]

    for raw_genome in raw_genomes:
        genomes.append([raw_genome[4], float(raw_genome[8]), float(raw_genome[11])])

    return genomes

def writeinfile(distance_matrix):
    dmfile = open("distancematrix.txt", "w")

    dmfile.write(distance_matrix.to_string())

    #for matrixline in distance_matrix[1:]:
    #    for linelement in matrixline[1:]:
    #        dmfile.write(str(linelement) + ' ')
    #    dmfile.write('\n\n')

    dmfile.close()
    return 0

genome_array = getgenomes()
genome_strings = [genome[0] for genome in genome_array]
distance_matrix = CreateDistanceMatrix(genome_strings)

writeinfile(distance_matrix)


