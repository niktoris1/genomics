from Get_Err_Rate import result, data, get_max_and_min_variants, stamms, start_time, end_time

outfile = open("result.txt", "w")

outfile.write('The sequencing error is: ' + str(result[1]) + '\n')

for read in data:
    outfile.write('Read of sample number ' + str(read.sample_id) + ' in position ' + str(read.position) + '\n')
    outfile.write(str(read.adenine_reads) + ' Adenine reads' + '\n')
    outfile.write(str(read.cytosine_reads) + ' Cytosine reads' + '\n')
    outfile.write(str(read.guanine_reads) + ' Guanine reads' + '\n')
    outfile.write(str(read.thymine_reads) + ' Thymine reads' + '\n')
    outfile.write('There are ' + str(read.number_of_variants) + ' valid variants' + '\n')

    if len(get_max_and_min_variants(read, read.number_of_variants)[0]):
        outfile.write(str(get_max_and_min_variants(read, read.number_of_variants)[0][0][0]) + ' is valid with share ' + str(read.share) + '\n\n')
    else:
        outfile.write(str(get_max_and_min_variants(read, read.number_of_variants)[0][0][0]) + ' is valid with share ' + str(read.share) + '\n\n')
        outfile.write(str(get_max_and_min_variants(read, read.number_of_variants)[0][1][0]) + ' is valid with share ' + str(1 - read.share) + '\n\n')

outfile.write('Overall results \n')

for sample in stamms:
    if sample[1] == sample[2]: # both variant are the same, which means that there is only one variant
        outfile.write('There is 1 stamm in sample ' + str(sample[0]) + '\n')
        outfile.write('The only one is ' + str(sample[1]) + '\n')
        outfile.write('The share is ' + str(sample[3]) + '\n\n')
    else:
        outfile.write('There are 2 stamms in sample ' + str(sample[0]) + '\n')
        outfile.write('First one is ' + str(sample[1]) + '\n')
        outfile.write('Second one is ' + str(sample[2]) + '\n')
        outfile.write('The share is ' + str(sample[3]) + '\n\n')


outfile.write('Time elapsed: ' + str(end_time - start_time) + ' seconds')
outfile.close()