#!/usr/bin/env python
"""
# Copyright 2015 Matt Settles
# Created Dec 1, 2014
"""
from optparse import OptionParser
import sys
import time
from subprocess import Popen, PIPE, STDOUT
import string
import re
import copy


def sp_gzip_read(file, bufsize=-1):
    p = Popen('gzip --decompress --to-stdout'.split() + [file], stdout=PIPE, stderr=STDOUT, bufsize=bufsize)
    return p.stdout


def sp_gzip_write(file, bufsize=-1):
    filep = open(file, 'wb')
    p = Popen('gzip', stdin=PIPE, stdout=filep, shell=True, bufsize=bufsize)
    return p.stdin

rcs = string.maketrans('TAGCtagc', 'ATCGATCG')


def revcomp(seq):
    return seq.translate(rcs)[::-1]


def rev(seq):
    return seq[::-1]


def Illumina_Sanger(qual):
    """
    Converts Illumina (Phred+64) to Sanger(Phred+33)
    """
    QualityScoreOut = ''
    for quality in qual:
        newQual = chr((ord(quality) - 64) + 33)
        QualityScoreOut += newQual
    return QualityScoreOut


def Illumina_ID(rid):
    """
    Converts BGI ids to Illumina 1.8 IDs
    """
    index = rid.find(":")  # finds the first occurance of ':'
    new_id = rid[:index] + ":1:12345" + rid[index:]
    new_id_split = re.split("#|/", new_id)
    new_id = new_id_split[0] + " " + new_id_split[2] + ":Y:0:" + new_id_split[1]
    return new_id


class fastqIter:
    " A simple file iterator that returns 4 lines for fast fastq iteration. "
    def __init__(self, handle):
        self.inf = handle

    def __iter__(self):
        return self

    def next(self):
        lines = {'id': self.inf.readline().strip(),
                 'seq': self.inf.readline().strip(),
                 '+': self.inf.readline().strip(),
                 'qual': self.inf.readline().strip()}
        assert(len(lines['seq']) == len(lines['qual']))
        if lines['id'] == '' or lines['seq'] == '' or lines['+'] == '' or lines['qual'] == '':
            raise StopIteration
        else:
            return lines

    @staticmethod
    def parse(handle):
        return fastqIter(handle)

    def close(self):
        self.inf.close()


def writeFastq(handle, fq):
    handle.write(fq['id'] + '\n')
    handle.write(fq['seq'] + '\n')
    handle.write(fq['+'] + '\n')
    handle.write(fq['qual'] + '\n')


def main(read1, read2, outfile1, outfile2, verbose):
    # Set up the global variables
    global pair_count
    global stime
    # Process Paired inputs:
    if read1.split(".")[-1] == "gz":
        iterator1 = fastqIter(sp_gzip_read(read1))
        # assume both are gz
        iterator2 = fastqIter(sp_gzip_read(read2))
    else:
        iterator1 = fastqIter(open(read1, 'r'))
        iterator2 = fastqIter(open(read2, 'r'))
    try:
        while 1:
            seq1 = iterator1.next()
            seq2 = iterator2.next()
            pair_count += 1
            # fix quality scores
            seq1['qual'] = Illumina_Sanger(seq1['qual'])
            seq2['qual'] = Illumina_Sanger(seq2['qual'])

            # @FCC76PFACXX:1:1101:1380:2035#ACTTGAATC_TCTTTCCCT/2
            seq1['id'] = Illumina_ID(seq1['id'])
            seq2['id'] = Illumina_ID(seq2['id'])

            writeFastq(outfile1, seq1)
            writeFastq(outfile2, seq2)
            if pair_count % 500000 == 0 and verbose:
                print "Pairs:", pair_count, "| reads/sec:", round(pair_count / (time.time() - stime), 0)
    except StopIteration:
        if verbose:
            print "Pairs:", pair_count, "| reads/sec:", round(pair_count / (time.time() - stime), 0)
        pass


#####################################
# Parse options and setup #
usage = "usage %prog -o [output file prefix (path + name)] --quite -1 [read1] -2 [read2]"
usage += "convertBGI2Illumina.py will process read files produced by BGI and make then more similar to current fastq files"
parser = OptionParser(usage=usage, version="%prog 0.0.1")

parser.add_option('-o', '--output', help="Directory + prefix to output reads",
                  action="store", type="str", dest="output_dir", default="paired-reads")

parser.add_option('-1', '--Read1', help="Read1, check length and ouput",
                  action="store", type="str", dest="read1", default=None)

parser.add_option('-2', '--Read2', help="Read2, check length and output",
                  action="store", type="str", dest="read2", default=None)
parser.add_option('-g', '--gzip', help="gzip the output",
                  action="store_true", dest="gzip", default=True)

parser.add_option('--quite', help="turn off verbose output",
                  action="store_false", dest="verbose", default=True)

(options, args) = parser.parse_args()

output_dir = options.output_dir

infile1 = options.read1
if infile1 is None:
    sys.stdout.write("Paired end file 1 is missing\n")
    sys.exit(1)
infile2 = options.read2
if infile2 is None:
    sys.stdout.write("Paired end file 2 is missing\n")
    sys.exit(1)
verbose = options.verbose


pair_count = 0
stime = time.time()

if options.gzip:
    outfile1 = sp_gzip_write(output_dir + "_R1_001.fastq.gz")
    outfile2 = sp_gzip_write(output_dir + "_R2_001.fastq.gz")
else:
    outfile1 = open(output_dir + "_R1_001.fastq", "w")
    outfile2 = open(output_dir + "_R2_001.fastq", "w")

main(infile1, infile2, outfile1, outfile2, verbose)

outfile1.close()
outfile2.close()

print "Final Pairs:", pair_count, "| reads/sec:", round(pair_count / (time.time() - stime), 0)

sys.exit(0)
