from collections import defaultdict
import gzip
import re
import subprocess
from itertools import groupby as g

old_illumina_header = ["instrument",
                       "flowcell_lane",
                       "flowcell_number",
                       None,  # x-tile
                       None,  # y-tile
                       "barcode",  # barcode - fetched later
                       "pair"]

illumina_header = ["instrument",
                   "run_id",
                   "flowcell_id",
                   "flowcell_lane",
                   None,  # tile number
                   None,  # x-tile
                   None,  # y-tile
                   "pair",
                   "filtered",
                   "control_bits",
                   "barcode"]  # barcode/index sequence; fetched later.

SRR_header = ['SRR']

stat_header = ["Total_Reads",
         "Unique_Reads",
         "Percent_Unique", 
         "Most_Abundant_Sequence",
         "Most_Abundant_Frequency",
         "Percentage_Unique_fq"]

def boolify(s):
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def autoconvert(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

def most_common(L):
    # Fetch most common item from a list.
  try:
    return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
  except:
    return ""

class fastq:
    # Simple class for reading fastq files.
    def __init__(self, filename):
        self.filename = filename
        # Get fastq information
        header_lines = [x["info"] for x in self.read(1)]
        header = re.split(r'(\:|#|/| )',header_lines[0])[::2]
        fetch_barcode = True

        if len(header) == 11:
            # Use new header format.
            use_header = illumina_header
        elif len(header) == 6:
            # Use old header format.
            use_header = old_illumina_header
        elif len(header) == 7:
            # Use old header and add pair if available.
            use_header = old_illumina_header + ["pair"]
        elif header[0].startswith("@SRR"):
            # Setup SRR Header
            header = header[0].split(".")[0:1]
            use_header = SRR_header
            fetch_barcode = False
        else:
            # If unknown header, enumerate
            use_header = ["h" + str(x) for x in range(0,len(header))]
            fetch_barcode = False

        if fetch_barcode == True:
            # Fetch index
            index_loc = use_header.index("barcode")
            fetch_index = [re.split(r'(\:|#|/| )',x["info"])[::2][index_loc] for x in self.read(1000)]
            self.barcode = most_common(fetch_index)


        self.header = {}
        # Set remaining attributes.
        for attr, val in zip(use_header, header):
            if attr is not None:
                val = autoconvert(val) # Set variable type
                self.header[attr] = val
                setattr(self, attr, val)

        # Fetch index
        #line["index"] = most_common(index_set)


    def read(self, n=-1):
        """
            Iterate through gzipped fastq file and put
            yield sequence+info in dictionary.
        """
        if self.filename.endswith(".gz"):
            open_file = gzip.open(self.filename, 'rb')
        else:
            open_file = open(self.filename, 'r') 
        with open_file as f:
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

    def calculate_fastq_stats(self):
        if self.filename.endswith(".fq"):
            # Read if not zipped.
            awk_read = "cat".format(**locals())
        else:
            # unzip
            awk_read = "gzcat"
        awk_one_liner =  """ awk '((NR-2)%4==0){ 
                                read=$1;total++;count[read]++
                             }
                             END{
                             for(read in count){
                             if(!max||count[read]>max) {
                                 max=count[read];
                                 maxRead=read};
                                 if(count[read]==1){
                                    unique++
                                 }
                             };
                             print total,
                                   unique,
                                   unique*100/total,
                                   maxRead,
                                   count[maxRead],
                                   count[maxRead]*100/total}'"""
        awk = awk_read + " " + self.filename + " | " + awk_one_liner
        stats = subprocess.check_output([awk], shell=True)
        stats = map(autoconvert,stats.strip().split(" "))
        stats = dict(zip(stat_header, stats))
        self.sequence_stats = stats
        return stats



x = fastq("../test/test1.fastq.gz")
x = fastq("../test/test2.fastq.gz")
x = fastq("../test/test3.fastq.gz")
x = fastq("../test/test3.fastq")
print x.header
