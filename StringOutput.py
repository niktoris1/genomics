from Get_Err_Rate import stamms

def GetStrings():

    class Genome:
        def __init__(self, genome_line, originating_sample, share_in_sample):
            self.genome_line = genome_line
            self.originating_sample = originating_sample
            self.share_in_sample = share_in_sample

    genomes = [] # genomes are organised as [line, sample_id, share]

    a = stamms

    for stamm in stamms:
        if stamm.dominant_stamm == stamm.non_dominant_stamm:
            genomes.append(Genome(stamm.dominant_stamm, stamm.sample_id, stamm.sample_share))
        else:
            genomes.append(Genome(stamm.dominant_stamm, stamm.sample_id, stamm.sample_share))
            genomes.append(Genome(stamm.non_dominant_stamm, stamm.sample_id, 1 - stamm.sample_share))

    outfile = open("resulting_genomes.txt", "w")

    for genome in genomes:
        outfile.write('There exist a genome ' + genome.genome_line + ' in a sample ' + str(genome.originating_sample) + ' with share ' + str(genome.share_in_sample) + '\n')

    outfile.close()

    return genomes

genomes = GetStrings()

