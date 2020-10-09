from Data_Read import data, get_max_and_min_variants, SNV_Reads
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
    true_variants = get_max_and_min_variants(read, 1)[0]  # returns dictionary with 1 name of the most frequent nycleotyde and number of reads with them
    #false_variants = get_max_and_min_variants(read, 1)[1] # returns dictionary with 3 names of other nycleotydes

    LLH_value = StirlingLogFactorial(read.total_coverage()) - StirlingLogFactorial(true_variants[0][1])- StirlingLogFactorial(read.total_coverage() - true_variants[0][1])+true_variants[0][1] * log(1 - error_rate)+(read.total_coverage() - true_variants[0][1]) * log(error_rate)

    return LLH_value


def LogLikelyhoodFunction2(read, error_rate, share): #returns LLH function with 2 true results. Shares is a share of the most frequent haplotype 0.5 < share < 1
    if error_rate <= 0 or error_rate >= 1:
        print('AWARE!!!')
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


def GetBestLLHValueByRead(read, error_rate, share): # We are given an error rate and share and have to give the answer - is there one or two maximums and if 2 (in assumption, that the share in 2 is known)? 
    LLH1_value = LogLikelyhoodFunction1(read, error_rate)
    LLH2_value = LogLikelyhoodFunction2(read, error_rate, share)

    if share >= 0.989: # We do not bother with samples, where second variant is extremely small
        return [1, LLH1_value]

    if LLH1_value > LLH2_value:
        return [1, LLH1_value] # return an indication, that LLH1 is better, plus LLH1 itself
    else:
        return [2, LLH2_value, share] # return an indication, that LLH2 is better, plus LLH2 itself and share, what share is the best
    
def ResultingLLHByPerson(sample_num, error_rate, share):
    LLH_Value = 0
    for read in data:
        if read.sample_id == sample_num:
            BestLLHValue = GetBestLLHValueByRead(read, error_rate, share)
            LLH_Value += BestLLHValue[1]
            read.number_of_variants = BestLLHValue[0]

            if read.number_of_variants == 2:
                read.share = BestLLHValue[2]
            else:
                read.share = 1
    return LLH_Value

def OptimiseLLHByPerson(sample_id, error_rate): # We optimise it in assumption, that we know an error

    start_share = 0.9
    LLH = scipy.optimize.minimize_scalar(
        lambda share, sample_id, error_rate: - ResultingLLHByPerson(sample_id, error_rate, share), bounds = (0, 1),
       args=(sample_id, error_rate), method='bounded', options={'xatol': 100000})

    LLH_value = - LLH.fun
    share = LLH.x

    print('LLH optimised for person', sample_id, 'value', LLH_value)
    
    return [LLH_value, share]

def ResultingLLHByData(data, error_rate): # get an LLH after optimisation by persons

    samples = []
    for read in data:
        if [read.sample_id] not in samples:
            samples.append([read.sample_id])

    LLH_value = 0
    for sample in samples:
        sample_res = OptimiseLLHByPerson(sample[0], error_rate)
        LLH_value += sample_res[0]
        sample.append(sample_res[1])

    return [LLH_value, samples] # return an array and samples in the following form [sample_id, share in sample]

def OptimiseLLHByData(data):

    start_error = 0.001
    LLH = scipy.optimize.minimize_scalar(
        lambda error_rate, data: - ResultingLLHByData(data, error_rate)[0], bounds=(0, 1),
        args=(data), method='bounded', options={'xatol': 100000})


    LLH_value = - LLH.fun
    LLH_error = LLH.x

    print('Value:', LLH_value, 'Error:', LLH_error)

    return [LLH_value, LLH_error]


def GetStamms(data):




    for sample in samples:
        sample.append('')
        sample.append('')

        sample_share = 1

        for read in data: # !!! in data we have more reads with som ids - why?
            if read.sample_id == sample[0]:
                if read.share == 1:
                    number_of_variants = 1
                else:
                    number_of_variants = 2

                if read.adenine_reads == 0 and read.cytosine_reads == 0 and read.guanine_reads == 0 and read.thymine_reads == 0:
                    sample[1] = sample[1] + '?'
                    sample[2] = sample[2] + '?'
                else:
                    variants = get_max_and_min_variants(read, number_of_variants)
                    if number_of_variants == 1:
                        sample[1] = sample[1] + GetLetter(variants[0][0][0])
                        sample[2] = sample[2] + GetLetter(variants[0][0][0])
                    if number_of_variants == 2:
                        sample[1] = sample[1] + GetLetter(variants[0][0][0])
                        sample[2] = sample[2] + GetLetter(variants[0][1][0])

                    if read.share != 1:
                        sample_share = read.share

        sample.append(sample_share)
        sample.append(len(sample[1]))


    return samples



start_time = time.time()

result = OptimiseLLHByData(data)
stamms = GetStamms(data)

end_time = time.time()








