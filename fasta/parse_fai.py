
def parse_fai(fai):
    """ Parses a fasta index file (produced by samtools index)
    and returns chromosome names and lengths """
    with open(fai) as f:
        return dict([x.strip().split("\t")[0:2] for x in f.readlines()])