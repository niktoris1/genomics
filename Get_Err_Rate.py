from Data_Read import data, get_max_and_min_variants, SNV_Reads
import math
from math import factorial, log
import scipy
from scipy.special import logsumexp
from scipy import optimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt
import numpy as np
import time



def perm(N, k): # have to write it from scratch, because in-built scipy.special.perm method is buggy in scipy 1.5.2 and returns inf. Since scipy 1.6.0 is supported only by older versions of python - we do this stuff
    return math.factorial(N) // (math.factorial(k) * (math.factorial(N-k)))

def LogLikelyhoodFunction1 (read, error_rate): #returns LLH function with one true result
    if error_rate <= 0 or error_rate >= 1:
        print('AWARE!!!')
        raise ValueError
    true_variants = get_max_and_min_variants(read, 1)[0]  # returns dictionary with 1 name of the most frequent nycleotyde and number of reads with them
    false_variants = get_max_and_min_variants(read, 1)[1] # returns dictionary with 3 names of other nycleotydes

    #print('ERR =', error_rate)
    LLH_value = log(perm(read.total_coverage(), true_variants[0][1])) + true_variants[0][1] * log(1 - error_rate) + \
                (read.total_coverage() - true_variants[0][1]) * log(error_rate)
    return LLH_value


def LogLikelyhoodFunction2(read, error_rate, share): #returns LLH function with 2 true results. Shares is an share of the most frequent haplotype 0.5 < share < 1
    if error_rate == 0 or error_rate == 1:
        raise ValueError
    true_variants = get_max_and_min_variants(read, 2)[0]  # returns dictionary with 2 names of the most frequent nycleotydes and number of reads with them
    false_variants = get_max_and_min_variants(read, 2)[1] # returns dictionary with other 2 nycleotydes

    assumed_total_for_1st_nycleotyde = math.floor(read.total_coverage() * share)
    assumed_total_for_2nd_nycleotyde = read.total_coverage() - assumed_total_for_1st_nycleotyde

    prob_list = [] # list of logarithms of probabilities of all variants
    # Here we add all probabilities, where we have the resulting distribution.
    # We consider all variants, how are results are distributed between true sequencing and errors in sequensing.

    for first_true in range(0, true_variants[0][1] + 1):
        for second_true in range(0, true_variants[1][1] + 1):
            first_false = true_variants[1][1] - second_true
            second_false = true_variants[0][1] - first_true
            first_err = assumed_total_for_1st_nycleotyde - first_true - first_false
            second_err = assumed_total_for_2nd_nycleotyde - second_true - second_false
            if first_err >= 0 and second_err >=0: # checking the correctness
                log1 = log(factorial(assumed_total_for_1st_nycleotyde)) - log(factorial(first_true)) - log(factorial(first_false)) - \
                       log(factorial(first_err)) + (first_true) * log(1 - error_rate) + \
                    (first_false) * log(error_rate / 3) + (first_err) * log(2 * error_rate / 3)

                log2 = log(factorial(assumed_total_for_2nd_nycleotyde)) \
                       - log(factorial(second_true)) - log(factorial(second_false)) - \
                       log(factorial(second_err)) + (second_true) * log(1 - error_rate) + \
                    (second_false) * log(error_rate / 3) + (second_err) * log(2 * error_rate / 3)

                prob_list.append(log1 + log2) # append the log of probability of considered distribution
    if len(prob_list) == 0:
        LLH_value = 0
    else:
        LLH_value = logsumexp(prob_list) # we use this trick to find a log sum of exponents
    return LLH_value


def GetBestLLHValue(read, error_rate): # We are given an error rate and have to give the answer - is there one or two maximums and if 2 - what is the share if 2? LLH for 1 is known, for 2 - we have to optimise it by share
    LLH1_value = LogLikelyhoodFunction1(read, error_rate)
    #print('LLH1', LLH1_value)
    start_share = [0.9]
    #bounds = Bounds([0.5], [0.99])
    #LLH2 = scipy.optimize.minimize(lambda share, read, error_rate: (-1) * LogLikelyhoodFunction2(read, error_rate, share), start_share, args=(read, error_rate), method='trust-constr', bounds = bounds)

    LLH2 = scipy.optimize.minimize(
        lambda share, read, error_rate: - LogLikelyhoodFunction2(read, error_rate, share), start_share,
        args=(read, error_rate), method='Nelder-Mead')

    #LLH2 = scipy.optimize.minimize_scalar(
    #    lambda share, read, error_rate: - LogLikelyhoodFunction2(read, error_rate, share), bounds = (0, 1),
    #    args=(read, error_rate), method='bounded')

    LLH2_value = - LLH2.fun
    LLH2_arg = LLH2.x

    #print('LLH2', LLH2_value)

    if LLH2_arg > 0.99: # We do not bother with samples, where second variant is extremely small
        #print('!!!!!!!!!!!!!!!!')
        #print("LLH with one variant. Value is", LLH1_value)
        return [1, LLH1_value]

    if LLH1_value > LLH2_value:
        #print("LLH with one variant. Value is", LLH1_value)
        return [1, LLH1_value] # return an indication, that LLH1 is better, plus LLH1 itself
    else:
        #print("LLH with two variants. Value is", LLH2_value, 'share', LLH2_arg)
        return [2, LLH2_value, LLH2_arg] # return an indication, that LLH2 is better, plus LLH2 itself and share, what share is the best


def ResultingLLH(data, error_rate):
    LLH = 0
    for read in data:
        LLH += GetBestLLHValue(read, error_rate)[1]
        read.number_of_variants = GetBestLLHValue(read, error_rate)[0]
        if read.number_of_variants == 2:
            read.share = GetBestLLHValue(read, error_rate)[2]
        else:
            read.share = 1
    print('Resulting LLH', LLH, 'with error rate', error_rate)
    return LLH


def GetErrorAndSplitReads(data): # returns error rate
    start_error = 0.002
    result = scipy.optimize.minimize(
        lambda error_rate, data: - ResultingLLH(data, error_rate), start_error,
        args=(data), method='Nelder-Mead', options= {'xatol': 1e-6}).x
    start_time = time.perf_counter()

    #result = scipy.optimize.minimize_scalar(
    #    lambda error_rate, data: - ResultingLLH(data, error_rate), bounds = (0, 1),
    #    args=(data), method='bounded', options= {'xatol': 1e-6}).x

    end_time = time.perf_counter()

    print('Error rate is', result, 'Time is', end_time-start_time)

    return [result, data]

[result, data] = GetErrorAndSplitReads(data)

print(result, data)


#x = np.arange(0.001, 0.01, 0.001)

#y = []
#for value in x:
#    y.append(ResultingLLH(data, value))

#plt.plot(x, y)
#plt.show()







