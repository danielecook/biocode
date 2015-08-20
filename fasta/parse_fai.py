
def parse_fai(fai, weight = False):
    """ Parses a fasta index file (produced by samtools index)
    and returns chromosome names and lengths

    Optionally generate weights
    """
    with open(fai) as f:
        chrom_length = [x.strip().split("\t")[0:2] for x in f.readlines()]
        chrom_length = [[x[0],int(x[1])] for x in chrom_length]
        if weight == True:
            genome_length = sum([x[1] for x in chrom_length])
            return dict([x[0],1.0*x[1]/genome_length] for x in chrom_length)
        else:
            return dict(chrom_length)