from Data_Read import data, max_and_min_variants, SNV_Reads, samples
from math_funcs import StirlingLogFactorial, GetLetter
import math
from math import log
import scipy
from scipy.special import logsumexp
from scipy import optimize

import time

def LogLikelyhoodFunction1 (read, error_rate): #returns LLH function with one true result
    if error_rate <= 0 or error_rate >= 1:
        print('AWARE!!!')
        raise ValueError
    true_variants = max_and_min_variants(read, 1).max_variants()  # returns dictionary with 1 name of the most frequent nycleotyde and number of reads with them
    #false_variants = get_max_and_min_variants(read, 1)[1] # returns dictionary with 3 names of other nycleotydes

    LLH_value = StirlingLogFactorial(read.total_coverage()) - StirlingLogFactorial(true_variants[0][1])- StirlingLogFactorial(read.total_coverage() - true_variants[0][1])+true_variants[0][1] * log(1 - error_rate)+(read.total_coverage() - true_variants[0][1]) * log(error_rate)

    return LLH_value


def LogLikelyhoodFunction2(read, error_rate, share): #returns LLH function with 2 true results. Shares is a share of the most frequent haplotype 0.5 < share < 1
    if error_rate <= 0 or error_rate >= 1:
        print('AWARE!!!')
        raise ValueError
    true_variants = max_and_min_variants(read, 2).max_variants()  # returns dictionary with 2 names of the most frequent nycleotydes and number of reads with them
    #false_variants = max_and_min_variants(read, 2).min_variants() # returns dictionary with other 2 nycleotydes

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
                log1 = StirlingLogFactorial(assumed_total_for_1st_nycleotyde) - StirlingLogFactorial(first_true) - StirlingLogFactorial(first_false) - \
                       StirlingLogFactorial(first_err) + (first_true) * log(1 - error_rate) + \
                       (first_false) * log(error_rate / 3) + (first_err) * log(2 * error_rate / 3)

                log2 = StirlingLogFactorial(assumed_total_for_2nd_nycleotyde) \
                       - StirlingLogFactorial(second_true) - StirlingLogFactorial(second_false) - \
                       StirlingLogFactorial(second_err) + (second_true) * log(1 - error_rate) + \
                       (second_false) * log(error_rate / 3) + (second_err) * log(2 * error_rate / 3)

                prob_list.append(log1 + log2) # append the log of probability of considered distribution
    if len(prob_list) == 0:
        LLH_value = 0
    else:
        LLH_value = logsumexp(prob_list) # we use this trick to find a log sum of exponents
    return LLH_value

def GetBestLLHValueByRead(read, error_rate, share): # We are given an error rate and share and have to give the answer - is there one or two maximums(in assumption, that the share in 2 is known)?


    LLH1_value = LogLikelyhoodFunction1(read, error_rate)

    if share >= 0.99: # We do not care for stamms with extra small shares
        read.number_of_variants = 1
        read.share = 1
        read.LLH_value = LLH1_value
        return LLH1_value

    LLH2_value = LogLikelyhoodFunction2(read, error_rate, share)

    if LLH1_value > LLH2_value:
        read.number_of_variants = 1
        read.LLH_value = LLH1_value
        read.share = 1
        return LLH1_value
    else:
        read.number_of_variants = 2
        read.share = share
        read.LLH_value = LLH2_value
        return LLH2_value



def ResultingLLHBySample(sample_id, error_rate, share):
    LLH_Sample_Value = 0
    for read in data:
        if read.sample_id == sample_id:
            LLH_Sample_Value += GetBestLLHValueByRead(read, error_rate, share)
    return LLH_Sample_Value


def ResultingLLHByData(data, error_rate, share_array): # share array is a share array for all samples

    LLH_value = 0

    for sample in samples:
        sample_res = ResultingLLHBySample(sample, error_rate, share_array[samples.index(sample)])
        LLH_value += sample_res

    return LLH_value  # return an array and samples in the following form [sample_id, share in sample]

def OptimiseLLHByData(data):

    start_error = 0.001
    start_share_array = [0.9] * len(samples)
    start_parameters = [start_error] + start_share_array

    resultwrapper = lambda parameters, data: - ResultingLLHByData(data, parameters[0], parameters[1:]) # parameters is an array of type [error_rate, share]

    print('Optimisation started')

    LLH = scipy.optimize.minimize(
        resultwrapper, start_parameters, args=(data), method='Nelder-Mead', tol=1e-2)

    print('Optimisation finished')

    LLH_value = - LLH.fun
    LLH_error = LLH.x[0]
    LLH_shares = LLH.x[1:]

    print('Value:', LLH_value, 'Error:', LLH_error, 'Shares:', LLH_shares)

    class LLH_data:
        def __init__(self, LLH_value, error_rate, share_array):
            self.LLH_value = LLH_value
            self.error_rate = error_rate
            self.share_array = share_array

    Optimised_LLH = LLH_data(LLH_value, LLH_error, LLH_shares)

    return Optimised_LLH




def GetStamms(data):

    class sample_with_stamms:
        def __init__(self, sample_id, dominant_stamm, non_dominant_stamm, sample_share):
            self.sample_id = sample_id
            self.dominant_stamm = dominant_stamm
            self.non_dominant_stamm = non_dominant_stamm
            self.sample_share = sample_share

    sample_array = []

    for sample in samples:

        sample = sample_with_stamms(sample_id = sample, dominant_stamm = '', non_dominant_stamm = '', sample_share = 1)

        sample_array.append(sample)

        sample_share = 1 # setting a default value of sample_share

        for read in data:
            if read.sample_id == sample.sample_id:
                if read.share == 1:
                    number_of_variants = 1
                else:
                    number_of_variants = 2

                if read.adenine_reads == 0 and read.cytosine_reads == 0 and read.guanine_reads == 0 and read.thymine_reads == 0:
                    sample.dominant_stamm = sample.dominant_stamm + '?'
                    sample.non_dominant_stamm = sample.non_dominant_stamm + '?'
                else:
                    variants = [max_and_min_variants(read, number_of_variants).max_variants(), max_and_min_variants(read, number_of_variants).min_variants()]
                    if number_of_variants == 1:
                        sample.dominant_stamm = sample.dominant_stamm + GetLetter(variants[0][0][0])
                        sample.non_dominant_stamm = sample.non_dominant_stamm + GetLetter(variants[0][0][0])
                    if number_of_variants == 2:
                        sample.dominant_stamm = sample.dominant_stamm + GetLetter(variants[0][0][0])
                        sample.non_dominant_stamm = sample.non_dominant_stamm + GetLetter(variants[0][1][0])

                    if read.share != 1:
                        sample_share = read.share

        sample.sample_share = sample_share
        #sample.append(len(sample[1])) - not sure what was here

    return sample_array # returns an array of samples with corresponding stamms and shares


start_time = time.time()

result = OptimiseLLHByData(data)
stamms = GetStamms(data)

end_time = time.time()

print('Optimisation took ', end_time - start_time, ' seconds')








