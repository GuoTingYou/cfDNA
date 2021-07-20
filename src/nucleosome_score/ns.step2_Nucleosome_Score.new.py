#!/home/guotingyou/anaconda3/bin/python
### Filename: step2_Nucleosome_Score.py ###
# -*- Author: GuoTingYou -*- #
# -*- Coding: utf-8 -*- #

import os
import sys
import time
import argparse

import re
import gzip
from functools import wraps
from itertools import takewhile
from collections import Counter, defaultdict, namedtuple

import numpy as np
from scipy.signal import fftconvolve

### I/O Interface and Usage ###

def argParse(N):
    script = sys.argv[0]
    usage = "Not enough arguments were given. Expected {1}\n{0}".format(
            " ".join(["Usage: python", script]), N)
    description = "{0:^81s}\n{1:^81s}\n{2:^81s}\n{3:^81s}".format(
                  "Filename: {}".format(script),
                  "-*- Author: GuoTingYou -*-",
                  "-*- Version: youyou.tortoise 0.0.0 -*-",
                  "-*- Coding: utf-8 -*-")
    parser = argparse.ArgumentParser(
             prog = script, usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    help_info = {
        "i": "Input file: read_starts file.",
        "I": "Input path: read_starts in the path will be inputs.",
        "l": "Input list: list of read_starts file.",
        "O": "Output: Absolute output path.",
        "o": "Output: Prefix of output name.",
        "m":("How to merge NS for multiple inputs."
             "\n`gene` : Merge samples of same ethnity. [Population level]"),
    }
    arguments = {
        "v": dict(action = "version", version = "youyou.tortoise 0.0.0"),
        "i": dict(type=str, help=help_info["i"], nargs="?", default=[], action='append'),
        "I": dict(type=str, help=help_info["I"], nargs="?", default='None'),
        "l": dict(type=open,help=help_info["l"], nargs="?", default=None),
        "O": dict(type=str, help=help_info["O"], nargs="?", default=os.getcwd()),
        "o": dict(type=str, help=help_info["o"], nargs="?", default='output'),
        "m": dict(type=str, help=help_info["m"], nargs="?", default='individual',
                  choices=['individual', 'gene', 'sample']),
    }
    parser.add_argument("-v", "--version",  **arguments["v"])
    parser.add_argument("-i", "--input",    **arguments["i"])
    parser.add_argument("-I", "--in_path",  **arguments["I"])
    parser.add_argument("-l", "--list",     **arguments["l"])
    parser.add_argument("-O", "--out_path", **arguments["O"])
    parser.add_argument("-o", "--out_name", **arguments["o"])
    parser.add_argument("-m", "--merge",    **arguments["m"])
    if len(sys.argv) <= N:
        parser.print_help()
        raise IndexError("Not enough arguments were given.")
    return parser.parse_args()
try:
    ioargs = argParse(1)
    i_file = ioargs.input    # read_starts
    i_list = ioargs.list     # list of samples
    i_path = ioargs.in_path  # input files path
    o_path = ioargs.out_path # must given
    o_name = ioargs.out_name + '.ns.gz'
    i_list = [line.strip() for line in i_list] if i_list is not None \
        else os.listdir(i_path) if os.path.exists(i_path) \
        else i_file if len(i_file) > 0 \
        else print('Null Input, Input is required.')
except Exception as err:
    raise Exception(err)
else:
    o_file = os.path.join(o_path, o_name)
    WINDOW = 'ToBeUpdate' # update from file
    SMOOTH = 'ToBeUpdate' # update from file
    MERGE = ioargs.merge
    if MERGE not in ('individual', 'sample', 'gene'):
        raise ValueError('Parameter "MERGE" should in ["individual", "sample", "gene"]')
    all_genes = []

### Class and Function ###

class Read(namedtuple('Read', 'gene, chrom, tss, strand, start, length')):
    def apply(self, **kwargs):
        kwargs = {name: value if not callable(value) \
                   else value(getattr(self, name))
                  for name, value in kwargs.items()}
        return self._replace(**kwargs)

def process_file(i_file):
    global WINDOW, SMOOTH, all_genes, MERGE
    result = defaultdict(list)

    def update_param(param, pattern, line):
        pattern = re.compile(pattern)
        match = pattern.search(line)
        if match is not None:
            param = int(match.group(1))
        else:
            param = param
        return param

    Gene = namedtuple('Gene', 'chrom, tss, name')

    with gzip.open(i_file, 'rt') as f:
        headers = (line for line in takewhile(lambda line: line.startswith('##'), f))
        content = (line.strip().split('\t') for line in f if not line.startswith('#'))

        for line in headers:
            WINDOW = update_param(WINDOW, r'window: (\d+)', line)
            SMOOTH = update_param(SMOOTH, r'smooth: (\d+)', line)

        for read in map(Read._make, content):
            gene = Gene(read.chrom, read.tss, read.gene)
            read = read.apply(tss=int, start=int)
            result[gene].append(read.start)
            if MERGE == 'gene' and gene not in all_genes:
                all_genes.append(gene)

    return {gene:Counter(result[gene]) for gene in result}
         # {gene: Counter()}

def merge_data(data:dict):
    global MERGE, all_genes
    result = {}

    # Just change data format
    if MERGE == 'individual':
        for ID, sample in data.items():
            for gene, counter in sample.items():
                result['::'.join([ID, *gene])] = counter
        return result
             # {sample_chrom_TSS:Counter()}

    # MERGE all genes within a sample
    elif MERGE == 'sample':
        for ID, sample in data.items():
            temp = Counter()
            for counter in sample.values():
                temp += counter
            result[ID] = temp
        return result
             # {sample:Counter()}

    # MERGE the same gene across samples
    elif MERGE == 'gene':
        for sample in data.values():
            for gene in all_genes:
                result[gene] = result.get(gene, Counter()) \
                             + sample.get(gene, Counter())
        return result
             # {gene:Counter()}

def Nucleosome_Score(counter):
    global WINDOW, SMOOTH
    n = int(SMOOTH/2) if SMOOTH % 2 == 0 else int((SMOOTH-1)/2)
    relativeRS = [counter.get(pos, 0)
                  for pos in range(-WINDOW-n-93, WINDOW+n+94)]
    kernel = np.r_[np.ones(20), np.ones(147) * 1j, np.ones(20)]
    ns_array = fftconvolve(relativeRS, kernel, 'valid')
    ns_array = (ns_array.real + 1) / (ns_array.imag + 1) * 7.35
    return '|'.join(np.round(ns_array, 4).astype(str))

def data_process(i_list):
    global MERGE
    sample_names = (
        os.path.basename(i_file).split('.')[0]
        for i_file in i_list
        if i_file.endswith('rss.gz')
    )
    data = {
        name: process_file(i_file)
        for name, i_file in zip(sample_names, i_list)
        if i_file.endswith('rss.gz')
    }
    return merge_data(data)

def output(data, o_file):
    stdoutSave = sys.stdout
    sys.stdout = gzip.open(o_file, 'wt')
    print('## window: {}'.format(WINDOW))
    print('## smooth: {}'.format(SMOOTH))
    print('# Merge_{}'.format(MERGE), 'Nucleosome Score', sep = '\t')
    for index, counter in data.items():
        if isinstance(index, tuple):
            index = '::'.join(index)
        NS = Nucleosome_Score(counter)
        print(index, NS, sep = '\t')
    sys.stdout.close()
    sys.stdout = stdoutSave

def main():
    # Data Process
    print('Data Process')
    data = data_process(i_list)

    # Calculate Nucleosome Score
    print('Calculate Nucleosome Score')
    output(data, o_file)

### Main ###

if __name__ == "__main__":

    print(os.path.realpath(__file__))
    print("\n{:-^81s}".format(
          " Task Start At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))

    main()

    print("{:-^81s}\n".format(
          " Task End At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))
    print("Practice_makes_perfect")
