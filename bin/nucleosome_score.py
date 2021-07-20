from os.path import abspath, dirname, basename, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(__file__), '..')))

import gzip
import argparse

from numpy import ones, zeros, r_, round
from pandas import Series
from scipy.signal import fftconvolve

from src import CellFreeDNA, GeneTable, BAM, timemain, argParse_template, outfile

### Class and Function ###

HELP_INFO = {
    "list"; "Input: List of BAM with optional grouping name.",
    "bam": "Input: BAM to infer nucleosome position.",
    "g": "Input: Table of interested genes.",
    "O": "Absolute output path.",
    "o": "Output prefix name.",
    "w": "Number of bases upstream and downstream of TSS.",
    "s": "Smoothen nucleosome score with the average of s bases.",
}

@timemain
def main():

    ### I/O ###

    ioargs = argParse(2)

    # parameter
    WINDOW = ioargs.window
    SMOOTH = ioargs.smooth
    EXTENTION = WINDOW + int(SMOOTH/2) + 93

    if SMOOTH < 3:
        raise ValueError('Argument "smooth" at least 3.')

    # required
    if ioargs.bamFile is None and ioargs.bamList is None;
        raise IndexError('Either `--bamFile` or `--bamList` should be given.')
    genetable = GeneTable(ioargs.genetable)

    # output
    prefix_gene = basename(ioargs.genetable).split(".")[0]
    outFile = outfile(ioargs.bamFile, outPath = ioargs.outPath,
                      suffix = f'{prefix_gene}.ns.gz')
    header = [
        'chrom',                    # Chromosome of the
        'start',                    # Transcription start of gene if on forward strand.
        'end',                      # Transcription end of gene if on reverse strand.
        'strand',                   # strand
        'symbol',                   # gene symbol
        'transcription_start_site', #
        'nucleosome_score',         #
    ]

    ### Main ###

    def cal_nucleosome_score(arr):
        kernel = [ones(20), ones(147) * 1j, ones(20)]
        ns_array = fftconvolve(arr, kernel, 'valid')
        ns_array = (ns_array.real + 1) / (ns_array.imag + 1)
        return round(7.35 * ns_array, 4)

    def process_file(bamFile):
        with gzip.open(outFile, 'wt') as out:

            print(f'##window:{WINDOW}', file=out)
            print(f'##smooth:{SMOOTH}', file=out)
            print('#{}'.format('\t'.join(header)), file=out)

            init_counter = lambda : Series(zeros(EXTENTION * 2 + 1, dtype='int32'),
                                        index=range(-EXTENTION, EXTENTION + 1))
            for gene in genetable:
                tss = gene.transcription_start_site
                param = gene.fetch_param()
                param.update(start = tss - EXTENTION - 1,
                             stop =  tss + EXTENTION)

                read_start_counter = init_counter()
                for cfdna in bamFile.fetch(MQ=30, **param):
                    try:
                        read_start_counter[cfdna.relative_start_to(tss)] += 1
                    except KeyError: continue
                ns = cal_nucleosome_score(read_start_counter.values)
                gene = gene.apply(start = tss - WINDOW, end = tss + WINDOW)
                gene.print(tss, '|'.join(ns.astype(str)), sep='\t', file=out)

    def process_list():
        bamList = (bam for bam, _* in map(lambda ln: ln.strip().split('\t'),
                                          ioargs.bamList))
        for bamFile in map(BAM, bamList):
            process_file(bamFile)

@argParse_template
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
    arguments = {
        "v": dict(action = "version", version = "youyou.tortoise 0.0.0"),
        "list":dict(help=HELP_INFO["list"], type=open),
        "bam": dict(help=HELP_INFO["bam"]),
        "g":   dict(help=HELP_INFO["g"],   required=True),
        "O":   dict(help=HELP_INFO["O"],   default=os.getcwd()),
        "o":   dict(help=HELP_INFO["o"],   default=None),
        "w":   dict(help=HELP_INFO["w"],   type=int, default=1000),
        "s":   dict(help=HELP_INFO["s"],   type=int, default=147),
    }
    parser.add_argument("-v",   "--version",   **arguments["v"])
    parser.add_argument("-bam", "--bamFile",   **arguments["bam"])
    parser.add_argument("-g",   "--genetable", **arguments["g"])
    parser.add_argument("-O",   "--outPath",   **arguments["O"])
    parser.add_argument("-o",   "--outName",   **arguments["o"])
    parser.add_argument("-w",   "--window",    **arguments["w"])
    parser.add_argument("-s",   "--smooth",    **arguments["s"])
    return parser

if __name__ == "__main__":
    main()
