import math

class SNV_Reads:
    def __init__(self, adenine_reads, cytosine_reads, guanine_reads, thymine_reads, position, sample_id, number_of_variants, share):
        self.adenine_reads = adenine_reads
        self.cytosine_reads = cytosine_reads
        self.guanine_reads = guanine_reads
        self.thymine_reads = thymine_reads
        self.position = position
        self.sample_id = sample_id
        self.number_of_variants = number_of_variants
        self.share = share

    def total_coverage(self): # returns total coverage
        return self.adenine_reads + self.cytosine_reads + self.guanine_reads + self.thymine_reads

    def count_freqs(self): # returns frequencies of all variants
        return [self.adenine_reads / self.total_coverage(), self.cytosine_reads / self.total_coverage(),
                self.guanine_reads / self.total_coverage(), self.thymine_reads / self.total_coverage()]

def get_max_and_min_variants(read, number_of_variants): # returns number_of_variants of nucleotides with a maximum number of reads
    dict = {'adenine_reads': read.adenine_reads,
            'cytosine_reads': read.cytosine_reads,
            'guanine_reads': read.guanine_reads,
            'thymine_reads': read.thymine_reads}

    dict = sorted(dict.items(), key=lambda item: item[1], reverse=True)

    return [list(dict)[:number_of_variants], list(dict)[number_of_variants:]]

file = open("data.txt", "r")

raw_data = [line.split() for line in file]

data = []
 # needed for quite unelegant solution of several similar reads on one position - python methos "index" returns the first occurence in the list. Fixed

for position_read in raw_data[1:]:
    position_num = 1
    for sample_reads in position_read[1:]:
        freqs = sample_reads.split('_')
        data.append(SNV_Reads(freqs[0], freqs[1], freqs[2], freqs[3], position_read[0], raw_data[0][position_num], 'Unknown', 'Unknown'))
        position_num = position_num + 1


samples = []

for read in data:
    if [read.sample_id] not in samples:
        samples.append([read.sample_id])
for sample in samples:
    counter = 0
    sample.append(counter)
    for read in data:
        if read.sample_id == sample[0]:
            sample[1] = sample[1] + 1




for read in data:
    if read.adenine_reads == 'NA' or read.cytosine_reads == 'NA' or read.guanine_reads == 'NA' or read.thymine_reads == 'NA':
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







print(len(data))









#number_of_reads = []

#for read in data:
#    number_of_reads.append(read.total_coverage())

#threshold = sorted(number_of_reads)[math.floor(len(number_of_reads) * 0.99)] #set a threshold for the snv's with small number of reads

#print(len(data))

#for read in data:
#    if read.total_coverage() < threshold:
#        data.remove(read) # remove all reads smaller than threshold

#print(len(data))

