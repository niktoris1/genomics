

class SNV_Reads:
    def __init__(self, adenine_reads, cytosine_reads, guanine_reads, thymine_reads, position, sample_id):
        self.adenine_reads = adenine_reads
        self.cytosine_reads = cytosine_reads
        self.guanine_reads = guanine_reads
        self.thymine_reads = thymine_reads
        self.position = position
        self.sample_id = sample_id



    def count_freqs(self):
        sum_reads = self.adenine_reads + self.cytosine_reads + self.guanine_reads + self.thymine_reads
        return [self.adenine_reads / sum_reads, self.cytosine_reads / sum_reads, self.guanine_reads / sum_reads, self.thymine_reads / sum_reads]

file = open("data.txt", "r")

raw_data = [line.split() for line in file]

data = []

for position_reads in raw_data[1:]:
    for sample_reads in position_reads[1:]:
        freqs = sample_reads.split('_')
        data.append(SNV_Reads(freqs[0], freqs[1], freqs[2], freqs[3], position_reads[0], raw_data[0][position_reads.index(sample_reads)]))

for read in data:
    if read.adenine_reads == 'NA' or read.cytosine_reads == 'NA' or read.guanine_reads == 'NA' or read.thymine_reads == 'NA':
        data.remove(read)

for read in data:
        read.adenine_reads = int(read.adenine_reads)
        read.cytosine_reads = int(read.cytosine_reads)
        read.guanine_reads = int(read.guanine_reads)
        read.thymine_reads = int(read.thymine_reads)
        read.position = int(read.position)
        read.sample_id = int(read.sample_id)

def percent_of_seldom_SNVs(read):
    true_read = max(read.adenine_reads, read.cytosine_reads, read.guanine_reads, read.thymine_reads)
    other_reads = read.adenine_reads + read.cytosine_reads + read.guanine_reads + read.thymine_reads - true_read
    return other_reads / (true_read + other_reads)

freqs = []

for read in data:
    freqs.append(percent_of_seldom_SNVs(read))

freqs.sort()


print(freqs)