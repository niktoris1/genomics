import numpy as np

logFactorialArray = [0] * 10000  # what should be instead of 10000

lf = 0
for k in range(1, len(logFactorialArray)):
    lf += np.log(k)
    logFactorialArray[k] = lf


def LogFactorial(k):  # Calculating log(n!)
    return logFactorialArray[k]

