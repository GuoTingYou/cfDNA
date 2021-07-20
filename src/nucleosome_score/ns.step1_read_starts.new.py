#!/home/guotingyou/anaconda3/bin/python
### Filename: step1_read_starts.py ###
# -*- Author: GuoTingYou -*- #
# -*- Coding: utf-8 -*- #

import os
import sys
import time

import re
import gzip
import argparse
from functools import wraps
from collections import namedtuple

import numpy as np
import pysam

### Class and Function ###

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
        "i": "Input file: BAM to pileup.",
        "g": "Input file: List of interest genes.",
        "O": "Absolute output path.",
        "o": "Output prefix name.",
        "w": "Number of bases upstream and downstream of TSS.",
        "s": "Smoothen nucleosome score with the average of s bases.",
    }
    arguments = {
        "v": dict(action = "version", version = "youyou.tortoise 0.0.0"),
        "i": dict(help=help_info["i"], required=True),
        "g": dict(help=help_info["g"], required=True),
        "O": dict(help=help_info["O"], default=os.getcwd()),
        "o": dict(help=help_info["o"], default=None),
        "w": dict(help=help_info["w"], type=int, default=1000),
        "s": dict(help=help_info["s"], type=int, default=147),
    }
    parser.add_argument("-v", "--version",   **arguments["v"])
    parser.add_argument("-i", "--input_bam", **arguments["i"])
    parser.add_argument("-g", "--gene_list", **arguments["g"])
    parser.add_argument("-O", "--out_path",  **arguments["O"])
    parser.add_argument("-o", "--out_name",  **arguments["o"])
    parser.add_argument("-w", "--window",    **arguments["w"])
    parser.add_argument("-s", "--smooth",    **arguments["s"])
    if len(sys.argv) <= N:
        parser.print_help()
        raise IndexError("Not enough arguments were given.")
    return parser.parse_args()

class BamRead:
    """
    Extend pysam read by delegate class.
    """
    def __init__(self, read):
        self.read = read
        self.flag = self.read.flag
        self.chrom = self.read.reference_name
        self.length = self.read.query_alignment_length
        self.strand = '-' if read.is_reverse else '+'
        self.start = self._get_read_start()

    @property
    def is_filter_pass(self, q=30):
        if (self.read.mapping_quality >= q and
        not self.read.is_supplementary and
        not self.read.is_duplicate and
        not self.read.is_secondary and
        not self.read.is_qcfail):
            return True
        return False

    @property
    def is_cigar_allM(self):
        # consider soft- and hard-clip
        pattern = re.compile(r'^\d+M$')
        cigar = self.read.cigarstring
        if cigar is None: return False
        return True if pattern.match(cigar) else False

    def start_relative_to(self, gene):
        """
        Adjust read start according to the gene's transcription start site.

        The result tell how many base pairs the read start site is upstream 
        or downstream relative to the TSS.
        """
        return self.start - gene.tss

    def _get_read_start(self):
        """
        Get read start site.
        If all bases of the read is aligned to reference (all Match/Mismatch),
        determined by CIGAR, start is the terminal of the read.
        Otherwise, start is the first aligned base of the read.

        Returns
        -------
        start : [int] read start site
            The position of terminal aligned base.
        """
        if self.is_cigar_allM:
            start = self.read.reference_end + 1 if self.read.is_reverse \
               else self.read.reference_start + 1 # 0 to 1-based
        else:
            align = self.read.query_alignment_end - 1 if self.read.is_reverse \
               else self.read.query_alignment_start # -1 for end point exculsive
            start = self._get_aligned_start(align)
        return start

    def _get_aligned_start(self, align):
        """
        The function is called when any base of the read is not aligned
        to the reference. It return the reference position of the last
        aligned base.

        Parameters
        ----------
        read : [pysam.Alignment object] (required)

        align : [int] first aligned position of the query (required)

        Returns
        -------
        rpos : [int, 1-based] reference position of the last aligned base
        """
        for qpos, rpos in self.read.get_aligned_pairs():
            if rpos and qpos == align:
               return rpos + 1 # to 1-based
        return 0 # something wrong

class GeneTable:
    Gene = namedtuple('Gene', 'gene, chrom, tss, pupS, pupE')

    def __init__(self, geneFile):
        self._file = self._open(geneFile)

    def __iter__(self):
        for gene in self._file:
            yield self._parse(gene)

    def _open(self, geneFile):
        return (line.strip().split('\t')
                for line in open(geneFile, 'r')
                if not line.startswith('#'))

    def _parse(self, gene):
        _, geneSymbol, chrom, start, end, strand = gene
        #chrom, start, end, strand, geneSymbol, *_ = gene
        chrom = self._parse_chrom(chrom)
        tss = self._get_tss(start, end, strand)
        pupS, pupE = self._get_pileup_region(tss)
        return self.Gene(geneSymbol, chrom, tss, pupS, pupE)

    def _parse_chrom(self, chrom):
        chrom = chrom if 'chr' in chrom else 'chr' + chrom
        return chrom.split('_')[0]

    def _get_tss(self, s, e, strand):
        return int(s) if strand in ('+', 'plus') \
          else int(e) if strand in ('-', 'minus') \
          else self.throw_error(f'Invalid {strand} for "strand"')

    def _get_pileup_region(self, tss):
        n = int(SMOOTH/2) if SMOOTH % 2 == 0 else int((SMOOTH-1)/2)
        pileup_start = tss - WINDOW - n - 93
        pileup_end   = tss + WINDOW + n + 93
        return pileup_start, pileup_end

    def throw_error(self, err: str):
        raise ValueError(err)

### Main ###

if __name__ == "__main__":

    print(os.path.realpath(__file__))
    print("\n{:-^81s}".format(
          " Task Start At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))

    ### I/O ###

    ioargs = argParse(2)
    # parameter
    WINDOW = ioargs.window
    SMOOTH = ioargs.smooth
    if SMOOTH < 3:
        raise ValueError('Argument "smooth" at least 3.')
    # required
    i_file = ioargs.input_bam
    i_gene = ioargs.gene_list
    # output
    o_path = ioargs.out_path
    o_name = ioargs.out_name if ioargs.out_name \
        else os.path.basename(i_file).split('.')[0]
    o_file = os.path.join(o_path, f'{o_name}.rss.gz')

    with pysam.AlignmentFile(i_file, 'rb') as bamFile, \
         gzip.open(o_file, 'wt') as o_file:

        header = '\t'.join('#Gene Chrom TSS Strand Start Length'.split(' '))
        print('## window: {}'.format(WINDOW), file=o_file)
        print('## smooth: {}'.format(SMOOTH), file=o_file)
        print(header, file=o_file)
        for gene in GeneTable(i_gene):
            #print(gene.gene, gene.chrom, gene.pupS, gene.pupE)
            #continue
            fetch_param = {
                'contig': gene.chrom,
                'start':  gene.pupS - 1, # pileup start
                'stop':   gene.pupE,     # pileip end
            }
            reads = (read for read in map(BamRead, bamFile.fetch(**fetch_param))
                          if read.is_filter_pass)
            for read in reads:
                relStart = read.start_relative_to(gene)
                print(gene.gene, read.chrom, gene.tss,
                      read.strand, relStart, read.length,
                      sep='\t', file=o_file)

    print("{:-^81s}\n".format(
          " Task End At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))
    print("Practice_makes_perfect")
