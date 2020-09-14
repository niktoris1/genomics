from Get_Err_Rate import result, data, get_max_and_min_variants

outfile = open("result.txt", "w")

outfile.write('The sequencing error is: ' + str(result) + '\n')

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

outfile.close()