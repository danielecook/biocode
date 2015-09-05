# biocode
Collection of scripts/utilities

__/fasta/fasta.py__ 

Utilities for working with reference genomes. Fasta class with methods for retrieving sequences from samtools indexed referenced genomes.


__/fasta/parse_fai.py__ - Parser for fasta index (__fai__) files.

__/fastq/fastq.py__ - A python fastq object. 

* Parse Header - Parses the first read and stores header information.
* __count_non_overlapping__ - Method for counting non-overlapping copies of a k-mer.
* __longest_sequential_repeat__ - Count frequency of longest sequential repeats.
* __calculate_fastq_stats__ - Aggregate statistics

