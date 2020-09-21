from math_funcs import StirlingLogFactorial

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 07:38:03 2020

@author: Dmitry
"""
import numpy as np
from numpy import random
from scipy import optimize

import time

random.seed(2039681)

# haps = [["A","A","T","C","A","G"],["A","A","G","C","A","G"],["G","C","G","T","T","A"]]
# haps = [[1,1,1,0,0,0],[1,1,0,0,0,0],[0,0,0,1,1,1]]
haplotypes = [[1, 1, 0, 0],
              [0, 0, 1, 1]]  # Ground truth haplotypes
haplotypeLength = len(haplotypes[0])
coverage = 100  # [80 + random.binomial(40, 0.5) for _ in hapLen] - coverage
num_of_haplotypes = len(haplotypes)  # Number of haplotypes
number_of_samples = 100  # Number of samples
er_rate = 0.01  # Error rate

class Sample:
    def __init__(self, num_of_haplotypes, coverage, haps, er_rate):
        self.prev = self.PrevSim(num_of_haplotypes)
        self.reads = self.SimReads(num_of_haplotypes, coverage, haps, er_rate)  # Observed data - read counts
        self.coef = []
        for snv in range(len(self.reads)):
            lf = StirlingLogFactorial(self.reads[snv][0] + self.reads[snv][1]) - StirlingLogFactorial(
                self.reads[snv][0]) - StirlingLogFactorial(self.reads[snv][1])
            self.coef.append(lf)
        self.coef = sum(self.coef)
        #print("the sum is", self.coef)

    # Simulate prevalence of each of num_of_haplotypes strains
    def PrevSim(self, num_of_haplotypes): #currently gives us a random prevalence of one strain - [0, 1] and [1, 0] with 0.3 chance and random distribition with 0.7 chance
        prevs = [0 for _ in range(num_of_haplotypes)]
        rnd = random.random()
        th = 0.9
        if num_of_haplotypes == 2:
            th = 1.0
        if rnd < 0.3:
            prevs[random.randint(0, 2)] = 1
        elif rnd < th:
            pops = random.choice([i for i in range(num_of_haplotypes)], 2, replace=False)
            prev = random.random()
            prevs[pops[0]] = prev
            prevs[pops[1]] = 1 - prev
        else:
            pops = random.choice([i for i in range(num_of_haplotypes)], 3, replace=False)
            prev1 = random.random()
            prev2 = prev1 * random.random()
            prevs[pops[0]] = prev1
            prevs[pops[1]] = prev2
            prevs[pops[2]] = 1 - prev1 - prev2
        return prevs

    # Simulate read counts for each SNV - single-nycleotid variant
    def SimReads(self, num_of_haplotypes, coverage, haps, er_rate):
        reads = []
        for snv in range(len(haps[0])):
            lcoverage = coverage - 20 + random.binomial(40, 0.5) # coverage is somewhat different from the real one
            read_counts = random.multinomial(lcoverage, pvals=self.prev)
            als = [0, 0]
            for hap_id in range(len(haps)):
                allele = haps[hap_id][snv]
                ers = random.binomial(read_counts[hap_id], er_rate)
                als[allele] += read_counts[hap_id] - ers
                als[1 - allele] += ers
            reads.append(als)
        # print('READS ARE', reads)
        return reads


# Ignore binomial coefficients, because they are constant through the whole run
# Notice, that it is MINUS log-llh
def sampleLLH(prev, haps, sample, er_rate):
    prev = np.append(prev, 1.0 - prev)
    llh = sample.coef
    if False:
        print(prev)
        print(haps)
        print(sample.reads)
    for snv in range(len(haps[0])):
        p = 0.0
        for i in range(len(prev)):
            p += prev[i] * ((1.0 - er_rate) * haps[i][snv] + er_rate * (1.0 - haps[i][snv]))
        if p <= 0:
            print(er_rate)
            print(prev)
            print([haps[i][snv] for i in range(len(prev))])
        if sample.reads[snv][0] != 0:
            llh += sample.reads[snv][0] * np.log(1.0 - p)

        if sample.reads[snv][1] != 0:
            llh += sample.reads[snv][1] * np.log(p)
    return -llh


def LLH(haps, num_of_haplotypes, er_rate, samples, prevs):
    haps = np.reshape(haps, (num_of_haplotypes, -1))
    print(haps)
    llh = 0

    for sample_num in range(len(samples)):
        llh += sampleLLH(prevs[sample_num], haps, samples[sample_num], er_rate)
        ++sample_num
    print('LLH:', llh)
    return llh

# Notice, that it is MINUS log-llh
def LLH_with_opt(haplotypes, num_of_haplotypes, er_rate, samples):
    haplotypes = np.reshape(haplotypes, (num_of_haplotypes, -1))
    print(haplotypes)
    llh = 0
    prevs = []
    for sample_num in range(len(samples)):
        haplotype_prevalence = [1 / num_of_haplotypes for _ in range(num_of_haplotypes - 1)]
        bounds = ((0, 1),) * (num_of_haplotypes - 1)
        #        bnds = optimize.Bounds([0, -0.5], [1.0, 2.0])
        linear_constraint = optimize.LinearConstraint([[1.0 for _ in range(num_of_haplotypes - 1)]], [0.0],
                                                      [1.0])  # for all x, 0 <= x <=1
        opt_res = optimize.minimize(sampleLLH, haplotype_prevalence, args=(haplotypes, samples[sample_num], er_rate), method='trust-constr',
                                    constraints=linear_constraint, bounds=bounds)

        prevs.append(opt_res.x)
        llh += opt_res.fun
    #print('LLH:', llh)
    return llh


def optimise_prevs_const_haps(haps, num_of_haplotypes, er_rate, samples, prevs):
    haps = np.reshape(haps, (num_of_haplotypes, -1))
    #print(haps)
    llh = 0

    out_prevs = []


    for sample_num in range(len(samples)):
        bnds = ((0, 1), ) * 1 # length of prev

        opt_res = optimize.minimize(sampleLLH, prevs[sample_num], args=(haps, samples[sample_num], er_rate), method='trust-constr', bounds=bnds, options = {'xtol': 1e-4, 'gtol': 1e-4, 'maxiter': 100})

        out_prevs.append(opt_res.x)
        llh += opt_res.fun

    #print('opt prevs:', out_prevs)
    #print('LLH:', LLH(haps, num_of_haplotypes, er_rate, samples, out_prevs))

    print ('OPTIMIZED PREVS')

    return out_prevs

def optimise_haps_const_prevs(haps, num_of_haplotypes, er_rate, samples, prevs):
    #haplotypes = np.reshape(haplotypes, (num_of_haplotypes, -1))

    bnds = ((0, 1), ) * len(haps)
    #linear_constraint = optimize.LinearConstraint([[1.0 for _ in range(len(bounds))]], [0.0],
    #                                                [1.0])      # for all x, 0 <= x <=1

    opt_res = optimize.minimize(LLH, haps, args=(num_of_haplotypes, er_rate, samples, prevs), method='trust-constr', bounds = bnds, options = {'xtol': 1e-4, 'gtol': 1e-4, 'maxiter': 100})
    opt_haps = opt_res.x

    #print('opt haps:', opt_haps)
    #print('LLH:', LLH(opt_haps, num_of_haplotypes, er_rate, samples, prevs))
    print ('OPTIMIZED HAPS')


    return opt_haps



def SimulateSamples(num_of_samples, num_of_haplotypes, coverage, haps, er_rate, sample=Sample):
    samples = [sample(num_of_haplotypes, coverage, haps, er_rate) for _ in range(num_of_samples)]
    return samples

samples = SimulateSamples(number_of_samples, num_of_haplotypes, coverage, haplotypes, er_rate)

start_time = time.time()

#haps_init = [0., 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
haps_init = [0.5 for _ in range(num_of_haplotypes * haplotypeLength)]
prev_init = [0.5 for _ in range(len(samples))]
haps_bnds = ((0, 1),) * len(haps_init)
prev_bnds = ((0, 1), )

opt_haps = optimise_haps_const_prevs(haps_init, num_of_haplotypes, er_rate, samples, prev_init)
opt_prevs = optimise_prevs_const_haps(opt_haps, num_of_haplotypes, er_rate, samples, prev_init)

for _ in range(5):
    opt_haps = optimise_haps_const_prevs(opt_haps, num_of_haplotypes, er_rate, samples, opt_prevs)
    opt_prevs = optimise_prevs_const_haps(opt_haps, num_of_haplotypes, er_rate, samples, opt_prevs)

print ("time elapsed: {:.2f}s". format(time.time() - start_time))

print ('RESULT')
print (opt_haps)
print (LLH(opt_haps, num_of_haplotypes, er_rate, samples, opt_prevs))





