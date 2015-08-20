from random import sample
from collections import OrderedDict
from bisect import bisect_left

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
            return OrderedDict([x[0],1.0*x[1]/genome_length] for x in chrom_length)
        else:
            return OrderedDict(chrom_length)

def generate_random_sites(fai, num_sites = 1000):
    """
        Generates random sites from a reference genome without replacement.
    """
    fai = parse_fai(fai)
    # Generate breakpoints
    breakpoints = []

    start = 0
    for k,v in fai.items():
        breakpoints.append(start + v)
        start = start + v
    genome_length = sum(fai.values())
    chromosomes = fai.keys()
    sites = sorted(sample(xrange(genome_length), num_sites))
    chromosomes =  [chromosomes[bisect_left(breakpoints, x)] for x in sites]
    # Process sites
    for k, site in enumerate(sites):
        bp = bisect_left(breakpoints, site)
        if bp > 0:
            sites[k] = site - breakpoints[bp-1]
    return zip(chromosomes, sites)
    