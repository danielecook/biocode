from collections import defaultdict
import gzip
import re
from itertools import groupby as g

old_illumina_header = [("instrument", str),
                       ("flowcell_lane", int),
                       ("flowcell_number", int),
                       None,  # x-tile
                       None,  # y-tile
                       ("barcode",str),  # barcode - fetched later
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
                   ("barcode",str)]  # barcode/index sequence; fetched later.

def most_common(L):
    # Fetch most common item from a list.
  try:
    return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
  except:
    return ""

def set_type(var, val):
    # Set variable types.
    if var == str:
        return str(val)
    elif var == int:
        return int(val)

class fastq:
    # Simple class for reading fastq files.
    def __init__(self, filename):
        self.filename = filename
        # Get fastq information
        header_lines = [x["info"] for x in self.read(1)]
        header = re.split(r'(\:|#|/| )',header_lines[0])[::2]
        if len(header) == 11:
            # Use new header format.
            use_header = illumina_header
        elif len(header) == 6:
            # Use old header format.
            use_header = old_illumina_header
        elif len(header) == 7:
            # Use old header and add pair if available.
            use_header = old_illumina_header + ["pair"]
        else:
            raise Exception("Unknown header")
        # Fetch index
        index_loc = use_header.index(("barcode",str))
        fetch_index = [re.split(r'(\:|#|/| )',x["info"])[::2][index_loc] for x in self.read(1000)]
        self.barcode = most_common(fetch_index)

        # Set remaining attributes.
        for attr, val in zip(use_header, header):
            if attr is not None:
                val = set_type(attr[1], val) # Set variable type
                setattr(self, attr[0], val)

        # Fetch index
        #line["index"] = most_common(index_set)


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
                if linenum < n or n == -1:
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
