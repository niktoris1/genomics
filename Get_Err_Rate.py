from Data_Read import data

number_of_reads = []

for read in data:
    number_of_reads.append(read.reads_total())

threshold = sorted(number_of_reads)[math.floor(len(number_of_reads) * 0.5)] #set a threshold for the snv's with small number of reads

for read in data:
    if read.reads_total() < threshold:
        data.remove(read) # remove all reads

