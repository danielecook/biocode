from collections import defaultdict
import gzip


class fastq:
    # Simple class for reading fastq files and performing simple operations.
    def __init__(self, filename):
        self.filename = filename

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
                if linenum <= n or n == -1:
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
