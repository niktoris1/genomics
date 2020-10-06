from Get_Err_Rate import result, data, get_max_and_min_variants, stamms, start_time, end_time



outfile = open("result.txt", "w")

outfile.write('The sequencing error is: ' + str(result[2]) + '\n')

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
    if len(sample) == 2:
        outfile.write('There is 1 stamm' + '\n')
        outfile.write('It is ' + str(sample[1]) + '\n\n')
    if len(sample) == 3:
        outfile.write('There are 2 stamms' + '\n')
        outfile.write('First one is ' + str(sample[1]) + '\n')
        outfile.write('Second one is ' + str(sample[2]) + '\n\n')

outfile.write('Time elapsed: ' + str(end_time - start_time) + ' seconds')
outfile.close()