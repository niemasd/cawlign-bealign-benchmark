#! /usr/bin/env python3
'''
Benchmark cawlign and bealign. Assumes both (and bam2msa) are in your PATH.
'''

# imports
from glob import glob
from gzip import open as gopen
from os import chdir, getcwd, mkdir
from os.path import abspath, expanduser, isfile, isdir
from random import sample
from shutil import rmtree
from subprocess import run
from sys import argv, stdout
from tqdm import tqdm
import argparse

# constants
DEFAULT_REPS = 10
OUTPUT_TMP_SUFFIX = '.tmp'
TIME_COMMAND_BASE = ['/usr/bin/time', '-v']
CAWLIGN_COMMAND_BASE = ['cawlign', '-r', 'HXB2_pol', '-t', 'codon', '-I', '-s', 'HIV_BETWEEN_F']
BEALIGN_COMMAND_BASE = ['bealign', '-r', 'HXB2_pol', '-m', 'HIV_BETWEEN_F']
BAM2MSA_COMMAND_BASE = ['bam2msa']
ZIP_COMMAND_BASE = ['zip', '-9', '-r']

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
    args.sequences = abspath(expanduser(args.sequences))
    assert isfile(args.sequences), "File not found: %s" % args.sequences
    args.output = abspath(expanduser(args.output))
    assert not isfile(args.output) and not isdir(args.output), "Output file exists: %s" % args.output
    args.output_tmp = '%s%s' % (args.output, OUTPUT_TMP_SUFFIX)
    assert not isdir(args.output_tmp) and not isfile(args.output_tmp), "Temporary output folder exists: %s" % args.output_tmp
    assert args.min_n > 0, "min_n (%s) must be positive" % args.min_n
    assert args.max_n > args.min_n, "max_n (%s) must be greater than min_n (%s)" % (args.max_n, args.min_n)
    assert args.delta > 0, "delta (%s) must be positive" % args.delta
    assert ((args.max_n - args.min_n) % args.delta) == 0, "max_n (%s) - min_n (%s) is not evenly divisible by delta (%s)" % (args.max_n, args.min_n, args.delta)
    return args

# read FASTA stream and return (ID,seq) dictionary
def readFASTA(fn):
    if fn.endswith('.gz'):
        stream = gopen(fn, 'rt')
    else:
        stream = open(fn, 'r')
    seqs = dict(); name = None; seq = ''
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
    IDs = sample(list(seqs.keys()), n)
    if out_fn.lower().endswith('.gz'):
        out_f = gopen(out_fn, 'wt')
    else:
        out_f = open(out_fn, 'w')
    for ID in IDs:
        out_f.write('>%s\n%s\n' % (ID, seqs[ID]))
    out_f.close()

# calculate Hamming distance (as a count)
def hamming_distance(s, t):
    assert len(s) == len(t), "Lengths of strings differ: %s and %s" % (len(s), len(t))
    return sum(1 for i in range(len(s)) if s[i] != t[i])

# main program
if __name__ == "__main__":
    # parse args and load sequences
    args = parse_args()
    seqs = readFASTA(args.sequences)
    assert args.max_n <= len(seqs), "max_n (%s) is larger than the total number of sequences in file: %s" % (args.max_n, args.sequences)
    print("Benchmark Command: %s" % ' '.join(argv)); stdout.flush()

    # create output tmp folder and benchmark datasets
    print("Creating benchmark datasets..."); stdout.flush()
    mkdir(args.output_tmp)
    len_max_n = len(str(args.max_n)); len_reps = len(str(args.reps))
    for n in range(args.min_n, args.max_n+1, args.delta):
        n_s = str(n).zfill(len_max_n)
        curr_n_folder = '%s/n%s' % (args.output_tmp, n_s)
        mkdir(curr_n_folder)
        for r in range(1, args.reps+1):
            r_s = str(r).zfill(len_reps)
            curr_rep_name = 'n%s.r%s' % (n_s, r_s)
            curr_rep_folder = '%s/%s' % (curr_n_folder, curr_rep_name)
            mkdir(curr_rep_folder)
            curr_rep_fasta_path = '%s/%s.fas' % (curr_rep_folder, curr_rep_name)
            subsample_seqs(seqs, n, curr_rep_fasta_path)

    # run bealign and cawlign on each dataset
    print("Running bealign and cawlign..."); stdout.flush()
    for seq_fn in tqdm(sorted(glob('%s/*/n*.r*/*.fas' % args.output_tmp))):
        seq_prefix = seq_fn.rstrip('.fas')

        # run cawlign
        cawlign_prefix = '%s.cawlign' % seq_prefix
        cawlign_time_fn = '%s.time.txt' % cawlign_prefix
        cawlign_aln_fn = '%s.aln' % cawlign_prefix
        cawlign_stdout_fn = '%s.out.txt' % cawlign_prefix
        cawlign_stderr_fn = '%s.err.txt' % cawlign_prefix
        cawlign_command = TIME_COMMAND_BASE + ['-o', cawlign_time_fn] + CAWLIGN_COMMAND_BASE + ['-o', cawlign_aln_fn, seq_fn]
        cawlign_stdout_f = open(cawlign_stdout_fn, 'w')
        cawlign_stderr_f = open(cawlign_stderr_fn, 'w')
        run(cawlign_command, stdout=cawlign_stdout_f, stderr=cawlign_stderr_f)
        cawlign_stdout_f.close()
        cawlign_stderr_f.close()

        # run bealign
        bealign_prefix = '%s.bealign' % seq_prefix
        bealign_time_fn = '%s.time.txt' % bealign_prefix
        bealign_bam_fn = '%s.bam' % bealign_prefix
        bealign_stdout_fn = '%s.out.txt' % bealign_prefix
        bealign_stderr_fn = '%s.err.txt' % bealign_prefix
        bealign_command = TIME_COMMAND_BASE + ['-o', bealign_time_fn] + BEALIGN_COMMAND_BASE + ['-K', seq_fn, bealign_bam_fn]
        bealign_stdout_f = open(bealign_stdout_fn, 'w')
        bealign_stderr_f = open(bealign_stderr_fn, 'w')
        run(bealign_command, stdout=bealign_stdout_f, stderr=bealign_stderr_f)
        bealign_stdout_f.close()
        bealign_stderr_f.close()

        # run bam2msa (convert bealign BAM output to FASTA)
        bam2msa_prefix = '%s.bam2msa' % bealign_prefix
        bam2msa_time_fn = '%s.time.txt' % bam2msa_prefix
        bam2msa_aln_fn = '%s.aln' % bam2msa_prefix
        bam2msa_stdout_fn = '%s.out.txt' % bam2msa_prefix
        bam2msa_stderr_fn = '%s.err.txt' % bam2msa_prefix
        bam2msa_command = TIME_COMMAND_BASE + ['-o', bam2msa_time_fn] + BAM2MSA_COMMAND_BASE + [bealign_bam_fn, bam2msa_aln_fn]
        bam2msa_stdout_f = open(bam2msa_stdout_fn, 'w')
        bam2msa_stderr_f = open(bam2msa_stderr_fn, 'w')
        run(bam2msa_command, stdout=bam2msa_stdout_f, stderr=bam2msa_stderr_f)
        bam2msa_stdout_f.close()
        bam2msa_stderr_f.close()

        # compare the alignments
        hamming_fn = '%s.hamming.tsv' % seq_prefix
        seqs_cawlign = readFASTA(cawlign_aln_fn)
        seqs_bealign = readFASTA(bam2msa_aln_fn)
        IDs_cawlign = sorted(seqs_cawlign.keys())
        IDs_bealign = sorted(seqs_bealign.keys())
        assert IDs_cawlign == IDs_bealign, "Mismatch between IDs in cawlign (%s) and bealign (%s)" % (cawlign_aln_fn, bam2msa_aln_fn)
        hamming_f = open(hamming_fn, 'w')
        hamming_f.write("ID\tHamming (count)\tLength\n")
        for ID in IDs_cawlign:
            s_cawlign = seqs_cawlign[ID]
            s_bealign = seqs_bealign[ID][:-3] # bealign includes terminal STOP in reference, so remove it
            d = hamming_distance(s_cawlign, s_bealign)
            hamming_f.write('%s\t%s\t%s\n' % (ID, d, len(s_cawlign)))
        hamming_f.close()

    # zip output directory and clean up
    print("Zipping benchmark results..."); stdout.flush()
    tmp = getcwd(); chdir(args.output_tmp)
    zip_command = ZIP_COMMAND_BASE + [args.output] + sorted(glob('*'))
    run(zip_command); chdir(tmp); rmtree(args.output_tmp)
