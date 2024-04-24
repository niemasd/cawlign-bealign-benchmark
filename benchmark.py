#! /usr/bin/env python3
'''
Benchmark cawlign and bealign. Assumes both are in your PATH.
'''

# imports
from gzip import open as gopen
from os.path import isfile
import argparse

# read FASTA stream and return (ID,seq) dictionary
def readFASTA(stream):
    seqs = {}
    name = None
    seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq
    return seqs

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sequences', required=True, type=str, help="Input Sequences (FASTA)")
    parser.add_argument('-n', '--min_n', required=True, type=int, help="Minimum Number of Sequences in Benchmark (inclusive)")
    parser.add_argument('-N', '--max_n', required=True, type=int, help="Maximum Number of Sequences in Benchmark (inclusive)")
    parser.add_argument('-d', '--delta', required=True, type=int, help="Delta Between Subsequent Valus of n")
    parser.add_argument('-r', '--reps', required=False, type=int, default=1, help="Number of Replicates per n")
    args = parser.parse_args()
    assert isfile(args.sequences), "File not found: %s" % args.sequences
    assert args.min_n > 0, "min_n (%s) must be positive" % args.min_n
    assert args.max_n > args.min_n, "max_n (%s) must be greater than min_n (%s)" % (args.max_n, args.min_n)
    assert args.delta > 0, "delta (%s) must be positive" % args.delta
    assert ((args.max_n - args.min_n) % args.delta) == 0, "max_n (%s) - min_n (%s) is not evenly divisible by delta (%s)" % (args.max_n, args.min_n, args.delta)
    return args

# main program
if __name__ == "__main__":
    # parse args and load sequences
    args = parse_args()
    if args.sequences.lower().endswith('.gz'):
        seqs_f = gopen(args.sequences, 'rt')
    else:
        seqs_f = open(args.sequences, 'r')
    seqs = readFASTA(seqs_f)
    assert args.max_n <= len(seqs), "max_n (%s) is larger than the total number of sequences in file: %s" % (args.max_n, args.sequences)
    pass # TODO subsample `seqs` with `reps` replicates for {min_n, min_n+delta, ..., max_n}