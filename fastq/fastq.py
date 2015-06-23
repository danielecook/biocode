from collections import defaultdict
import gzip
import re

old_illumina_header = [("instrument", str),
                       ("flowcell_lane", str),
                       ("flowcell_number", int),
                       None,  # x-tile
                       None,  # y-tile
                       ("index", str),
                       ("pair", int)]

illumina_header = [("instrument", str),
                   ("run_id", int),
                   ("flowcell_id", str),
                   ("flowcell_lane", int),
                   None,  # tile number
                   None,  # x-tile
                   None,  # y-tile
                   ("pair", int),
                   ("filtered", str),
                   ("control_bits", int),
                   None]  # index sequence

def set_type(var, val):
    if var == str:
        return str(val)
    elif var == int:
        return int(val)

class fastq:
    # Simple class for reading fastq files.
    def __init__(self, filename):
        self.filename = filename
        # Get fastq information
        header_lines = [x["info"] for x in self.read(100)]
        header = re.split(r'(\:|#| )',header_lines[0])[::2]
        if len(header) == 11:
            self.instrument = header[0]
            for attr, val in zip(illumina_header, header):
                if attr is not None:
                    val = set_type(attr[1], val) # Set variable type
                    setattr(self, attr[0], val)


        elif len(header) == 21:
            pass
        else:
            raise Exception("Unknown header")

    def read(self, n=-1):
        """
            Iterate through gzipped fastq file and put
            yield sequence+info in dictionary.
        """
        with gzip.open(self.filename, 'rb') as f:
            dna = {}
            for linenum, line in enumerate(f):
                dna["info"] = line.strip()
                dna["seq"] = f.next().strip()
                f.next()
                dna["qual"] = f.next().strip()
                if linenum*4 <= n or n == -1:
                    yield dna
                else:
                    break

    def count_non_overlapping(self, needle, n=-1):
        """
            Count the number of non-overlapping repeats
            and return a tuple of the distribution.
        """
        count = defaultdict(int)
        for i in self.read(n):
            count[i["seq"].count(needle)] += 1
        return sorted(count.items())

    def longest_sequential_repeat(self, needle, n=-1):
        """
            Count the longest sequential repeat,
            and return a tuple of the distribution.
        """
        count = defaultdict(int)
        for i in self.read(n):
            needle_rep = needle
            while True:
                if i["seq"].find(needle_rep) > 0:
                    needle_rep += needle
                else:
                    break
            count[len(needle_rep)/len(needle)-1] += 1
        return sorted(count.items())
