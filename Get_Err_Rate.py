from Data_Read import data, SNV_Reads, get_max_and_min_variants
import math

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
    print (true_variants, false_variants)

    LLH_value = 0

    good_variant = 0
    for variant in true_variants:
        good_variant += variant[1]

    LLH_value += math.log(perm(read.total_coverage(), good_variant)) + good_variant * (1 - error_rate) + (read.total_coverage() - good_variant) * (error_rate)

    print(LLH_value)


def LogLikelyhoodFunction2 (read, error_rate, share): #returns LLH function with 2 true results. Shares is an share of the most frequent haplotype 0.5 < share < 1
    if error_rate == 0 or error_rate == 1:
        raise ValueError
    true_variants = get_max_and_min_variants(read, 2)[0]  # returns dictionary with 2 names of the most frequent nycleotydes and number of reads with them
    false_variants = get_max_and_min_variants(read, 2)[1] # returns dictionary with other 2 nycleotydes
    print (true_variants, false_variants)

    LH_value = 0


    assumed_total_for_1st_nycleotyde = math.floor(read.total_coverage() * share)
    assumed_total_for_2nd_nycleotyde = read.total_coverage() - assumed_total_for_1st_nycleotyde

    for first_true in range(0, true_variants[0][1] + 1):
        second_false = true_variants[0][1] - first_true
        for second_true in range(0, true_variants[1][1] + 1):
            first_false = true_variants[1][1] - second_true
            if first_true + first_false <= assumed_total_for_1st_nycleotyde and second_true + second_false <= assumed_total_for_2nd_nycleotyde: # checking the correctness
                LH_value += perm(assumed_total_for_1st_nycleotyde, first_true) * (1 - error_rate) ** first_true * (error_rate) ** (assumed_total_for_1st_nycleotyde - first_true) + \
                            perm(assumed_total_for_1st_nycleotyde - first_true, first_false) * (0.5) ** first_false + \
                            perm(assumed_total_for_2nd_nycleotyde, second_true) * (1 - error_rate) ** second_true * (error_rate) ** (assumed_total_for_2nd_nycleotyde - second_true) + \
                            perm(assumed_total_for_2nd_nycleotyde - second_true, second_false) * (0.5) ** second_false

                a = perm(assumed_total_for_1st_nycleotyde, first_true) * (1 - error_rate) ** first_true * (error_rate) ** (assumed_total_for_1st_nycleotyde - first_true)
                b = perm(assumed_total_for_1st_nycleotyde - first_true, first_false) * (0.5) ** first_false
                c = perm(assumed_total_for_2nd_nycleotyde, second_true) * (1 - error_rate) ** second_true * (error_rate) ** (assumed_total_for_2nd_nycleotyde - second_true)
                d = perm(assumed_total_for_2nd_nycleotyde - second_true, second_false) * (0.5) ** second_false


                print('First true =', first_true, 'of', true_variants[0][1], '\n',
                    'Second false =', true_variants[0][1] - first_true, '\n',
                    'Second true =', second_true, 'of', true_variants[1][1], '\n',
                    'First false =', true_variants[1][1] - second_true, '\n')

    print(LH_value)


#for read in data:
#   LogLikelyhoodFunction2(read, 0.1, 0.9)

LogLikelyhoodFunction2(data[0], 0.1, 0.7)


