import math

class SNV_Reads: # class describes a sum of all reads in a certain position in a certain sample
    def __init__(self, adenine_reads, cytosine_reads, guanine_reads, thymine_reads, position, sample_id, number_of_variants, share, LLH_value):
        self.adenine_reads = adenine_reads
        self.cytosine_reads = cytosine_reads
        self.guanine_reads = guanine_reads
        self.thymine_reads = thymine_reads
        self.position = position # position of said read in a genome
        self.sample_id = sample_id # sample from which the reads are taken
        self.number_of_variants = number_of_variants # number of valid stamms in a read - 1 or 2
        self.share = share # share of dominant stamm. Value of 1 means that there is only 1 stamm there
        self.LLH_value = LLH_value

    def total_coverage(self): # returns total coverage
        return self.adenine_reads + self.cytosine_reads + self.guanine_reads + self.thymine_reads

    def count_freqs(self): # returns frequencies of all variants as an array
        return [self.adenine_reads / self.total_coverage(), self.cytosine_reads / self.total_coverage(),
                self.guanine_reads / self.total_coverage(), self.thymine_reads / self.total_coverage()]

class max_and_min_variants:
    def __init__(self, read, number_of_variants):
        self.read = read
        self.number_of_variants = number_of_variants

    def get_dict(self):
        dict = {'adenine_reads': self.read.adenine_reads,
            'cytosine_reads': self.read.cytosine_reads,
            'guanine_reads': self.read.guanine_reads,
            'thymine_reads': self.read.thymine_reads}

        sorted_dict = sorted(dict.items(), key=lambda item: item[1], reverse=True) # redefining dict

        return sorted_dict

    def max_variants(self):
        return list(dict)[:self.number_of_variants]

    def min_variants(self):
        return list(dict)[self.number_of_variants:]

def data_read():
    file = open("data.txt", "r")
    raw_data = [line.split() for line in file]
    data = []

    for position_read in raw_data[1:]:
        position_num = 1 # needed for quite unelegant solution of several similar reads on one position - python methos "index" returns the first occurence in the list. Fixed for now
        for sample_reads in position_read[1:]:
            freqs = sample_reads.split('_')
            data.append(SNV_Reads(freqs[0], freqs[1], freqs[2], freqs[3], position_read[0], raw_data[0][position_num], 'Unknown', 'Unknown', 'Undefined'))
            position_num = position_num + 1

    for read in data:
        if read.adenine_reads == 'NA' or read.cytosine_reads == 'NA' or read.guanine_reads == 'NA' or read.thymine_reads == 'NA': # if we have an uncertainity, we return zeros
            read.adenine_reads = 0
            read.cytosine_reads = 0
            read.guanine_reads = 0
            read.thymine_reads = 0

    for read in data:
            read.adenine_reads = int(read.adenine_reads)
            read.cytosine_reads = int(read.cytosine_reads)
            read.guanine_reads = int(read.guanine_reads)
            read.thymine_reads = int(read.thymine_reads)
            read.position = int(read.position)
            read.sample_id = int(read.sample_id)

    return data

def get_samples():

    samples = []
    for read in data:
        if [read.sample_id] not in samples:
            samples.append([read.sample_id])

    return samples

data = data_read()
samples = get_samples()



