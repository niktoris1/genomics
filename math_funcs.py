import math
from math import log

def StirlingLogFactorial(n):
    if n == 0:
        return 0
    else:
        return (n + 0.5) * log(n) - n + 0.5 * log(2 * math.pi)