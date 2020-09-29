import math
from math import log

def StirlingLogFactorial(n):
    if n == 0:
        return 0
    else:
        return (n + 0.5) * log(n) - n + 0.5 * log(2 * math.pi)

def GetLetter(variant):
    if variant == 'adenine_reads':
        answer = 'A'
    if variant == 'cytosine_reads':
        answer = 'C'
    if variant == 'guanine_reads':
        answer = 'G'
    if variant == 'thymine_reads':
        answer = 'T'

    return answer