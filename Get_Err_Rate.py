from Data_Read import data, get_max_and_min_variants
import math
from math import factorial, log
import scipy
from scipy.special import logsumexp
from scipy import optimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt
import numpy as np


number_of_reads = []

for read in data:
    number_of_reads.append(read.total_coverage())

threshold = sorted(number_of_reads)[math.floor(len(number_of_reads) * 0.5)] #set a threshold for the snv's with small number of reads

for read in data:
    if read.total_coverage() < threshold:
        data.remove(read) # remove all reads smaller than threshold

def perm(N, k): # have to write it from scratch, because in-built scipy.special.perm method is buggy in scipy 1.5.2 and returns inf. Since scipy 1.6.0 is supported only by older versions of python - we do this stuff
    return math.factorial(N) // (math.factorial(k) * (math.factorial(N-k)))

def LogLikelyhoodFunction1 (read, error_rate): #returns LLH function with one true result
    if error_rate == 0 or error_rate == 1:
        raise ValueError
    true_variants = get_max_and_min_variants(read, 1)[0]  # returns dictionary with 1 name of the most frequent nycleotyde and number of reads with them
    false_variants = get_max_and_min_variants(read, 1)[1] # returns dictionary with 3 names of other nycleotydes

    LLH_value = 0

    good_variant = 0
    for variant in true_variants:
        good_variant += variant[1]

    LLH_value += log(perm(read.total_coverage(), good_variant)) + good_variant * (1 - error_rate) + (read.total_coverage() - good_variant) * (error_rate)


    return LLH_value


def LogLikelyhoodFunction2(read, error_rate, share): #returns LLH function with 2 true results. Shares is an share of the most frequent haplotype 0.5 < share < 1
    if error_rate == 0 or error_rate == 1:
        raise ValueError
    true_variants = get_max_and_min_variants(read, 2)[0]  # returns dictionary with 2 names of the most frequent nycleotydes and number of reads with them
    false_variants = get_max_and_min_variants(read, 2)[1] # returns dictionary with other 2 nycleotydes

    LLH_value = 0

    assumed_total_for_1st_nycleotyde = math.floor(read.total_coverage() * share)
    assumed_total_for_2nd_nycleotyde = read.total_coverage() - assumed_total_for_1st_nycleotyde

    for first_true in range(0, true_variants[0][1] + 1): # here we add all probabilities, where we have the resulting distribution. We consider all variants, how are results are distributed between true sequencing and errors in sequensing.
        second_false = true_variants[0][1] - first_true
        for second_true in range(0, true_variants[1][1] + 1):
            first_false = true_variants[1][1] - second_true
            if first_true + first_false <= assumed_total_for_1st_nycleotyde and second_true + second_false <= assumed_total_for_2nd_nycleotyde: # checking the correctness
                log1 = log(factorial(assumed_total_for_1st_nycleotyde)) - log(factorial(first_true)) - log(factorial(first_false)) - \
                       log(factorial(assumed_total_for_1st_nycleotyde - first_true - first_false)) + first_true * (1 - error_rate) + \
                    (assumed_total_for_1st_nycleotyde - first_true) * (error_rate)

                log2 = log(factorial(assumed_total_for_2nd_nycleotyde)) - log(factorial(second_true)) - log(factorial(second_false)) - \
                       log(factorial(assumed_total_for_2nd_nycleotyde - second_true - second_false)) + second_true * (1 - error_rate) + \
                    (assumed_total_for_2nd_nycleotyde - second_true) * (error_rate)

                LLH_value = logsumexp([log1, log2]) # we use this trick to find a log og sum of exponents
    return LLH_value


def GetMinLLHValue(read, error_rate): # We are given an error rate and have to give the answer - is there one or two maximums and if 2 - what is the share? LLH for 1 is known, for 2 - we have to optimise it by share
    LLH1_value = LogLikelyhoodFunction1(read, error_rate)
    start_share = [0.75]
    #bounds = Bounds([0.5], [0.99])
    #LLH2 = scipy.optimize.minimize(lambda share, read, error_rate: (-1) * LogLikelyhoodFunction2(read, error_rate, share), start_share, args=(read, error_rate), method='trust-constr', bounds = bounds)
    LLH2 = scipy.optimize.minimize(
        lambda share, read, error_rate: - LogLikelyhoodFunction2(read, error_rate, share), start_share,
        args=(read, error_rate), method='Nelder-Mead')

    LLH2_value = - LLH2.fun
    LLH2_arg = LLH2.x

    if LLH2_arg > 0.99: # We do not bother with very small samples of second type
        print("LLH in assumption of one true variant:", LLH1_value)
        return [1, LLH1_value]

    if LLH1_value > LLH2_value:
        print("LLH in assumption of one true variant:", LLH1_value)
        return [1, LLH1_value] # return an indication, that LLH1 is better, plus LLH1 itself
    else:
        print("LLH in assumption of two true variants:", LLH2_value, 'share', LLH2_arg)
        return [2, LLH2_value, LLH2_arg] # return an indication, that LLH2 is better, plus LLH2 itself and share, where it is the best


def ResultingLLH(data, error_rate):
    LLH = 0
    for read in data:
        LLH += GetMinLLHValue(read, error_rate)[1]
        #read.number_of_variants = GetMinLLHValue(read, error_rate)[0]
        #if read.number_of_variants == 2:
        #    read.share = GetMinLLHValue(read, error_rate)[2]
        #else:
        #    read.share = 1
    print(LLH, error_rate)
    return LLH


def GetErrorAndSplitReads(data): # returns error rate
    start_error = 0.1
    result = scipy.optimize.minimize(
        lambda error_rate, data: - ResultingLLH(data, error_rate), start_error,
        args=(data), method='Nelder-Mead')

    print('Error rate is', result.x)

    return result.x


x = np.arange(0.01, 0.1, 0.01)

y = []
for value in x:
    y.append(ResultingLLH(data, value))

plt.plot(x, y)
plt.show()




