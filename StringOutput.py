from Get_Err_Rate import stamms

def GetStrings():

    genomes = [] # genomes are organised as [line, sample_id, share]

    for stamm in stamms:
        if stamm[1] == stamm[2]:
            genomes.append([stamm[1], stamm[0], stamm[3], len(stamm[1])])
        else:
            genomes.append([stamm[1], stamm[0], stamm[3], len(stamm[1])])
            genomes.append([stamm[2], stamm[0], 1 - stamm[3], len(stamm[1])])

    outfile = open("resulting_genomes.txt", "w")

    for genome in genomes:
        outfile.write('There exist a genome ' + genome[0] + ' in a sample ' + str(genome[1]) + ' with share ' + str(genome[2]) + '\n')

    return 0

GetStrings()
