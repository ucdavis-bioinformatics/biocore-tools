#!/usr/bin/env python

"""
Scripts to Analyzer a sam/bam file for quality reporting information not commonly found in samtools or bamtools

# Copyright 2016, Matt Settles
# Modified June 8, 2016
"""
import os
import sys
import numpy
import re
import warnings

from subprocess import Popen
from subprocess import PIPE
import signal
import shlex

'''
https://samtools.github.io/hts-specs/SAMv1.pdf

NOTES, SAM/BAM format
1 QNAME String
2 FLAG Int
3 RNAME String
4 POS Int
5 MAPQ Int
6 CIGAR String
7 RNEXT String
8 PNEXT Int
9 TLEN Int
10 SEQ String
11 QUAL String

Note (the & operation produces a value equivalent to the binary with a 1 in that position):
0x1 template having multiple segments in sequencing
0x2 each segment properly aligned according to the aligner
0x4 segment unmapped
0x8 next segment in the template unmapped
0x10 SEQ being reverse complemented
0x20 SEQ of the next segment in the template being reversed
0x40 the first segment in the template
0x80 the last segment in the template
0x100 secondary alignment
0x200 not passing quality controls
0x400 PCR or optical duplicate

## contig info in header
@SQ Reference sequence dictionary. The order of @SQ lines defines the alignment sorting order.
     SN* Reference sequence name. Each @SQ line must have a unique SN tag. The value of this field is used in the alignment records in RNAME and RNEXT fields. Regular expression: [!-)+-<>-~][!-~]*
     LN* Reference sequence length. Range: [1,231-1]
     AS Genome assembly identifier.
     M5 MD5 checksum of the sequence in the uppercase, excluding spaces but including pads (as '*'s).
     SP Species.
     UR URI of the sequence. This value may start with one of the standard protocols, e.g http: or ftp:. If it does not start with one of these protocols, it is assumed to be a file-system path.
'''


# Open up a bam file, process first with samtools
def sp_bam2sam(file):
    p = Popen(shlex.split('samtools view -h') + [file],
              stdout=PIPE,
              stderr=None,
              bufsize=-1,
              preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    if p.returncode:
        raise
    return p.stdout


# Process a Sam/Bam Cigar String
class cigarString:
    '''
    Process a sam file CIGAR string
    '''
    pattern = re.compile('([MIDNSHPX=])')

    def __init__(self, cigar):
        values = self.pattern.split(cigar)[:-1]
        self.paired = (values[n:n + 2] for n in xrange(0, len(values), 2))  # pair values by twos

    def getAlignmentLength(self):
        g = 0
        for pair in self.paired:
            l = int(pair[0])
            t = pair[1]
            if t == 'M':
                g += l
            elif t == 'I':
                pass
            elif t == 'D':
                g += l
            elif t == 'N':
                pass
            elif t == 'S':
                pass
            elif t == 'H':
                pass
            elif t == 'P':
                pass
            elif t == 'X':
                pass
            elif t == '=':
                pass
            else:
                warnings.warn("encountered unhandled CIGAR character %s" % t)
                pass
        return g


# Store Contig relavant data from SamBam
class contigData:
    '''
    class to store data associated with a contig for reporting purposes
    '''
    secondary_alignment_count = 0
    single_end_count = 0
    paired_end_count = 0
    paired_end_alignment_count = 0
    paired_end_broken = 0

    def __init__(self, contig_name, contig_length):
        self.contig_name = contig_name
        self.contig_length = contig_length

    def add_secondary_align(self):
        self.secondary_alignment_count += 1

    def add_single_count(self):
        self.single_end_count += 1

    def add_paired_count(self):
        self.paired_end_count += 1

    def add_paired_broken_count(self):
        self.paired_end_broken += 1

    def add_paired_end_align_count(self):
        self.paired_end_alignment_count += 1


# Primary processing function
def processSamBam(insam, output):
    ''' given a sam (or bam) file in insam (already open),
        process the file line by line generating relevant mapping statisitics
    '''
    # variables to store
    i = 0
    single_end_count = 0
    paired_end_count = 0

    # when encountering a PE1 or PE2 read, store until the matching read is found
    PE1 = {}
    PE2 = {}

    contig_data = {}

    for line in insam:
        # Comment/header lines start with @
        if i % 1000000 == 0 and i > 0:
            print "Records processed: %s, single end reads: %s, paired end reads: %s" % (i, single_end_count, paired_end_count)

        # Header LINES
        if line[0] == "@":  # header line
            if line[1:3] == "SQ":
                # initiate the contig
                line = line.strip().split('\t')[1:]
                tags = {}
                for pair in line:
                    pair = pair.split(':')
                    tags[pair[0]] = pair[1]
                try:
                    contig_data[tags['SN']] = contigData(tags['SN'], tags['LN'])
                except KeyError:
                    sys.stderr("header contig info did not contain required fields SN and/or LN")
                    sys.exit(1)
            else:  # skip header lines that are not contig info
                continue
        # Aligned read LINES
        elif len(line.strip().split()) >= 11:  # not header lines and contains minimum number of columns
            i += 1
            line = line.strip().split()
            contig = str(line[2])
            flag = int(line[1])

            if (flag & 0x100):  # secondary alignment
                contig_data[contig].add_secondary_align()
                continue
            # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1
            # Handle SE:
            if not (flag & 0x1) and not (flag & 0x4):  # If a single end read present, ignore it
                contig_data[contig].add_single_count()
                single_end_count += 1
                continue
            # Handle PE:
            # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
            if (flag & 0x1) and ((flag & 0x4) or (flag & 0x8)):
                # one of the two reads in the pair is unmapped
                if (flag & 0x40):
                    contig_data[contig].add_paired_count()
                    paired_end_count += 1
                    contig_data[contig].add_paired_broken_count()
                if (flag & 0x4):
                    contig_data[contig].add_paired_end_align_count()
                continue
            # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
            if (flag & 0x1) and not (flag & 0x4) and not (flag & 0x8):
                contig_data[contig].add_paired_end_align_count()
                ID = line[0].split("#")[0]
                cigar = cigarString(line[5])
                if (flag & 0x40):  # If read 1
                    contig_data[contig].add_paired_count()
                    paired_end_count += 1
                    if (flag & 0x10):  # if RevComp
                        r1 = ['R', line[2], int(line[3]) + cigar.getAlignmentLength() - 1, line[3]]
                    else:
                        r1 = ['F', line[2], line[3], int(line[3]) + cigar.getAlignmentLength() - 1]
                    if ID in PE2:
                        # write_pair_tab(r1, PE2[ID])
                        del PE2[ID]
                    else:
                        PE1[ID] = r1
                elif (flag & 0x80):  # If read 2 (last segment in template)
                    if (flag & 0x10):  # if RevComp
                        r2 = ['R', line[2], int(line[3]) + cigar.getAlignmentLength() - 1, line[3]]
                    else:
                        r2 = ['F', line[2], line[3], int(line[3]) + cigar.getAlignmentLength() - 1]
                    if ID in PE1:
                        # write_pair_tab(PE1[ID], r2)
                        del PE1[ID]
                    else:
                        PE2[ID] = r2
    # end sam/bam loop
    print "Records processed: %s, single end reads: %s, paired end reads: %s" % (i, single_end_count, paired_end_count)

# Main function
if len(sys.argv) == 2:
    infile = sys.argv[1]
    if not os.path.exists(infile):
        print "Error, can't find input file %s" % infile
        sys.exit(1)

    if infile.split(".")[-1] == "bam":
        insam = sp_bam2sam(infile)
    elif infile.split(".")[-1] == "sam":
        insam = open(infile, 'r')
    else:
        sys.stderr("Error, requires a sam/bam (ends in either .sam or .bam) file as input")
        sys.exit(1)
    filen, ext = os.path.splitext(infile)
    output = filen + "_report"
elif len(sys.argv) == 1:
    # reading/writing from stdin/stdout
    insam = sys.stdin
    output = "analyzeSamBam_report"
else:
    print "Usage: analyzeSamBam.py samfile.sam OR analyzeSamBam.py samfile.bam"
    sys.exit()


processSamBam(insam, output)

insam.close()




# if len(sys.argv) != 4:
#     print "Usage: produce_histograms.py file1.fastq.gz file2.fastq.gz file3_SE.fastq.gz"
#     sys.exit()

# fq1 = sys.argv[1]
# fq2 = sys.argv[2]
# fqse = sys.argv[3]

# # fq1 = "BS_Eggs_R1.fastq"
# # fq2 = "BS_Eggs_R2.fastq"
# # fqse = "BS_Eggs_SE.fastq"

# if fq1.split('.')[-1] == 'gz':
#     h1 = SeqIO.parse(gzip.open(fq1, 'rb'), 'fastq')
# else:
#     h1 = SeqIO.parse(open(fq1, 'r'), 'fastq')

# if fq2.split('.')[-1] == 'gz':
#     h2 = SeqIO.parse(gzip.open(fq2, 'rb'), 'fastq')
# else:
#     h2 = SeqIO.parse(open(fq2, 'r'), 'fastq')

# if fqse.split('.')[-1] == 'gz':
#     s1 = SeqIO.parse(gzip.open(fqse, 'rb'), 'fastq')
# else:
#     s1 = SeqIO.parse(open(fqse, 'r'), 'fastq')


# pairs_i = 0
# pairs_v = []
# try:
#     print "Analyzing Pairs"
#     while True:
#         r1 = h1.next()
#         r2 = h2.next()
#         pairs_v.extend([len(r1.seq) + len(r2.seq)])
#         pairs_i += 1
# except StopIteration:
#     pass
# print "Finished processing %s pairs of reads." % pairs_i 

# singles_i = 0
# singles_v = []
# try:
#     print "Analyzing Singles"
#     while True:
#         se1 = s1.next()
#         singles_v.extend([len(se1.seq)])
#         singles_i += 1
# except StopIteration:
#     pass
# print "Finished processing %s single reads." % singles_i 

# h1.close()
# h2.close()
# s1.close()

# print "Producing histogram files"
# puni = numpy.unique(pairs_v)
# phist, bin_edges = numpy.histogram(pairs_v, bins=numpy.append(puni,puni[-1]+1) - 0.1)

# suni = numpy.unique(singles_v)
# shist, bin_edges = numpy.histogram(singles_v, bins=numpy.append(suni,suni[-1]+1) - 0.1)

# def write_hist_file(hist_file, hist, unique):
#     fp = open(hist_file, 'w')
#     for i in range(min(unique), max(unique)+1):
#         try:
#             value = hist[numpy.where(unique==i)][0]
#         except (ValueError, IndexError):
#             value = ""
#         fp.write( "%s\t%s\n" % (i, value))
#     fp.close()


# def write_histogram_file(histogram_file, hist, unique):
#     max_num_asterisks = float(72)
#     scale = float(max_num_asterisks) / float(max(hist))
#     fp = open(histogram_file, 'w')
#     fp.write("# each asterisk represents approximately %s reads\n" % int(max(hist)/max_num_asterisks))
#     for i in range(min(unique), max(unique)+1):
#         try:
#             num_asterisks = int(scale*hist[numpy.where(unique==i)][0])
#         except (ValueError, IndexError):
#             num_asterisks = 0
#         fp.write( "%s\t%s\n" % (i, '*' * num_asterisks))
#     fp.close()

# write_hist_file(os.path.join(os.path.dirname(os.path.abspath(fq1)),"pairs.hist"), phist, puni)

# write_histogram_file(os.path.join(os.path.dirname(os.path.abspath(fq1)),"pairs.histogram"), phist, puni)

# write_hist_file(os.path.join(os.path.dirname(os.path.abspath(fqse)),"singles.hist"), shist, suni)

# write_histogram_file(os.path.join(os.path.dirname(os.path.abspath(fqse)),"singles.histogram"), shist, suni)
