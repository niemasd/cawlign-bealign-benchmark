#! /usr/bin/env python3
'''
Benchmark cawlign and bealign. Assumes both are in your PATH.
'''

# imports
from glob import glob
from gzip import open as gopen
from os import mkdir
from os.path import isfile, isdir
from random import choices
from subprocess import run
import argparse

# constants
DEFAULT_REPS = 10
OUTPUT_TMP_SUFFIX = '.tmp'
TIME_COMMAND_BASE = ['/usr/bin/time', '-v']
CAWLIGN_COMMAND_BASE = ['cawlign', '-r', 'HXB2_pol', '-t', 'codon', '-I', '-s', 'HIV_BETWEEN_F']

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sequences', required=True, type=str, help="Input Sequences (FASTA)")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output File (ZIP)")
    parser.add_argument('-n', '--min_n', required=True, type=int, help="Minimum Number of Sequences in Benchmark (inclusive)")
    parser.add_argument('-N', '--max_n', required=True, type=int, help="Maximum Number of Sequences in Benchmark (inclusive)")
    parser.add_argument('-d', '--delta', required=True, type=int, help="Delta Between Subsequent Valus of n")
    parser.add_argument('-r', '--reps', required=False, type=int, default=DEFAULT_REPS, help="Number of Replicates per n")
    args = parser.parse_args()
    assert isfile(args.sequences), "File not found: %s" % args.sequences
    assert not isfile(args.output) and not isdir(args.output), "Output file exists: %s" % args.output
    args.output_tmp = '%s%s' % (args.output, OUTPUT_TMP_SUFFIX)
    assert not isdir(args.output_tmp) and not isfile(args.output_tmp), "Temporary output folder exists: %s" % args.output_tmp
    assert args.min_n > 0, "min_n (%s) must be positive" % args.min_n
    assert args.max_n > args.min_n, "max_n (%s) must be greater than min_n (%s)" % (args.max_n, args.min_n)
    assert args.delta > 0, "delta (%s) must be positive" % args.delta
    assert ((args.max_n - args.min_n) % args.delta) == 0, "max_n (%s) - min_n (%s) is not evenly divisible by delta (%s)" % (args.max_n, args.min_n, args.delta)
    return args

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

# create subsampled sequence dataset
def subsample_seqs(seqs, n, out_fn):
    assert n < len(seqs), "n (%s) is larger than the number of sequences (%s)" % (n, len(seqs))
    IDs = choices(list(seqs.keys()), k=n)
    if out_fn.lower().endswith('.gz'):
        out_f = gopen(out_fn, 'wt')
    else:
        out_f = open(out_fn, 'w')
    for ID in IDs:
        out_f.write('>%s\n%s\n' % (ID, seqs[ID]))
    out_f.close()

# main program
if __name__ == "__main__":
    # parse args and load sequences
    args = parse_args()
    if args.sequences.lower().endswith('.gz'):
        seqs_f = gopen(args.sequences, 'rt')
    else:
        seqs_f = open(args.sequences, 'r')
    seqs = readFASTA(seqs_f); seqs_f.close()
    assert args.max_n <= len(seqs), "max_n (%s) is larger than the total number of sequences in file: %s" % (args.max_n, args.sequences)

    # create output tmp folder and benchmark datasets
    mkdir(args.output_tmp)
    len_max_n = len(str(args.max_n)); len_reps = len(str(args.reps))
    for n in range(args.min_n, args.max_n+1, args.delta):
        n_s = str(n).zfill(len_max_n)
        for r in range(1, args.reps+1):
            r_s = str(r).zfill(len_reps)
            curr_rep_name = 'n%s.r%s' % (n_s, r_s)
            curr_rep_folder = '%s/%s' % (args.output_tmp, curr_rep_name)
            mkdir(curr_rep_folder)
            curr_rep_fasta_path = '%s/%s.fas' % (curr_rep_folder, curr_rep_name)
            subsample_seqs(seqs, n, curr_rep_fasta_path)

    # run bealign and cawlign on each dataset
    for seq_fn in sorted(glob('%s/n*.r*/*.fas' % args.output_tmp)):
        seq_prefix = seq_fn.rstrip('.fas')

        # run cawlign
        cawlign_prefix = '%s.cawlign' % seq_prefix
        cawlign_command = TIME_COMMAND_BASE + ['-o', '%s.time.txt' % cawlign_prefix] + CAWLIGN_COMMAND_BASE + ['-o', '%s.aln' % cawlign_prefix, seq_fn]
        cawlign_stdout_f = open('%s.cawlign.out.txt' % seq_prefix, 'w')
        cawlign_stderr_f = open('%s.cawlign.err.txt' % seq_prefix, 'w')
        cawlign_stdout_f.write('COMMAND: %s\n\n' % ' '.join(cawlign_command))
        run(cawlign_command, stdout=cawlign_stdout_f, stderr=cawlign_stderr_f)
        cawlign_stdout_f.close(); cawlign_stderr_f.close()