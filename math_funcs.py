import math
from math import log

import numpy as np

def StirlingLogFactorial(n):
    if n == 0:
        return 0
    else:
        return (n + 0.5) * log(n) - n + 0.5 * log(2 * math.pi)

def GetLetter(variant):
    if variant == 'adenine_reads':
        answer = 'A'
    elif variant == 'cytosine_reads':
        answer = 'C'
    elif variant == 'guanine_reads':
        answer = 'G'
    elif variant == 'thymine_reads':
        answer = 'T'

    return answer

def HammingDistance(first_string, second_string):
    distance = 0
    for letter_num in range(len(first_string)):
        if first_string[letter_num] is not second_string[letter_num]:
            distance = distance + 1
    return distance

def CreateDistanceMatrix(string_array):
    distance_matrix = np.array([[0] * len(string_array)] * len(string_array))
    for first_string_num in range(len(string_array)):
        for second_string_num in range(len(string_array)):
            distance_matrix[first_string_num][second_string_num] = HammingDistance(string_array[first_string_num], string_array[second_string_num])
    return distance_matrix
