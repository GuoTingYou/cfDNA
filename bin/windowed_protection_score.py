from os.path import realpath, abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

import argparse
import bisect
import gzip
import pysam
from collections import defaultdict

from src import (
    # class
    BAM,
    GeneTable,
    # function
    argParse_template,
    timemain,
    outfile,
)
HELP_INFO = {
    "bam": "Input: <cfDNA.bam> with absolute path.",
    "bed": "Input:",
    "g":   "Input: <tissue_specific_gene> with absolute path.",
    "r":   "Chromosome to compare VCF.",
    "O":   "Output: Absolute path for output file.",
    "q":   "Parameter: Mapping Quality of READ in BAM above `q` will be pileuped.",
    "Q":   "Parameter: Base Quality of READ in BAM above `Q` will be pileuped.",
    "window_size": "Parameter:",
    "gene_format": "Parameter:",
    "exclude_supplementary": "Parameter:",
    "exclude_duplicate": "Parameter:",
    "exclude_secondary": "Parameter:",
    "proper_pair_only": "Parameter:",
}
@timemain
def main():

    ### I/O ###

    ioargs = argParse(1)

    # input
    CHROM = ioargs.chromosome
    bamFile = BAM(ioargs.bamFile)
    bedFile = ioargs.bedFile
    #genetable = GeneTable(ioargs.genetable, fmt=ioargs.gene_format)

    # output
    outFile = outfile(
        ioargs.bamFile, outPath = ioargs.outPath,
        suffix = f'wps.gz'
    )
    global HEADER, WPS_DICT
    HEADER = [
        'chrom',        # Chromosome of the SNP (with chr prefix).
        'start',        # Start position of the SNP (0-based coordinate).
        'end',          # End position of the SNP (end exclusive).
        'genotype',     # Combination genotype of the SNP locus.
        'allele_A',     # Share allele of combination genotype of the SNP.
        'allele_B',     # Special allele of combination genotype of the SNP.
        'fwdDepthA',    # Depth of forward oriented cfDNA with share allele.
        'rvsDepthA',    # Depth of reverse oriented cfDNA with share allele.
        'fwdDepthB',    # Depth of forward oriented cfDNA with special allele.
        'rvsDepthB',    # Depth of reverse oriented cfDNA with special allele.
        'fraction',     # Fraction of cfDNA originated from interested sample.
        'insertsize_A', # Insertsize of cfDNA with share allele.
        'insertsize_B', # Insertsize of cfDNA with special allele.
    ]
    WPS_DICT = defaultdict(lambda: defaultdict(lambda: 0))

    # parameters
    global WINDOW_SIZE
    filters = {
        'MQ': ioargs.mapq, # maping quality
        'exclude_supplementary': ioargs.exclude_supplementary, # default False
        'exclude_duplicate': ioargs.exclude_duplicate,         # default False
        'exclude_secondary': ioargs.exclude_secondary,         # default False
        'proper_pair_only': ioargs.proper_pair_only,           # default False
        'unique_score': 10,
        'mismatch': 5,
    }
    WINDOW_SIZE = ioargs.window_size

    ### main ###

    chromosomes = tuple(f'chr{N}' for N in [*range(1, 23), 'X'])

    with gzip.open(outFile, 'wt') as out:

        if bedFile is not None:
            bed = (l.strip().split('\t') for l in gzip.open(bedFile, 'rt'))
            for chrom, start, end, *_ in bed:
                fragments = get_fragments(
                    bamFile, chrom, int(start), int(end), **filters)
                renew_wps_dict(fragments, out)

        elif CHROM == 'genome':
            fragments = get_fragments(bamFile, **filters)
            renew_wps_dict(fragments, out)

        else:
            fragments = get_fragments(bamFile, CHROM, **filters)
            renew_wps_dict(fragments, out)

def get_fragments(bamFile, *args, **kwargs):
    f = (read.fragment() for read in bamFile.fetch(*args, **kwargs))
    return (fragment.typed for fragment in f if not fragment.is_sentinel)

def renew_wps_dict(fragments, out):
    global WPS_DICT
    for fragment in fragments:
        for pos, score in zip(
            fragment.windowed_protection_coord(window_size=WINDOW_SIZE),
            fragment.windowed_protection_array(window_size=WINDOW_SIZE)):
            WPS_DICT[fragment.chrom][pos] += score
        else:
            write(fragment, out)
    else:
        write(fragment, out, vacant=True)

def write(fragment, out, vacant=False):
    global WPS_DICT
    if vacant:
        for chrom in WPS_DICT:
            for pos in WPS_DICT[chrom]:
                wps = WPS_DICT[chrom].pop(pos)
                print(chrom, pos, pos+1, wps, sep='\t', file=out)

    for chrom in WPS_DICT:
        if chrom <= fragment.chrom:
            coords = list(WPS_DICT[chrom].keys())
            bound = int(fragment.start - WINDOW_SIZE / 2)
            while coords and coords[0] < bound:
                pos = coords.pop(0)
                wps = WPS_DICT[chrom].pop(pos)
                print(chrom, pos, pos+1, wps, sep='\t', file=out)

@argParse_template
def argParse(N):
    usage = (f"python {sys.argv[0]} "
              "-bam <cfDNA.bam> -bg <background.VCF> -it <interested.VCF> [-O ...]")
    description = f"Filename: {sys.argv[0]}"
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    arguments = {
        "v":   dict(action="version", version="youyou 0.0.0"),
        "bam": dict(help=HELP_INFO["bam"]),
        "bed": dict(help=HELP_INFO["bed"], default=None),
        "g":   dict(help=HELP_INFO["g"],   default=None),
        "r":   dict(help=HELP_INFO["r"], default='genome'),
        "O":   dict(help=HELP_INFO["O"], default=os.getcwd()),
        "q":   dict(help=HELP_INFO["q"], type=int, default=30),
        "Q":   dict(help=HELP_INFO["Q"], type=int, default=20),
        "window_size": dict(help=HELP_INFO["window_size"], type=int, default=120),
        "gene_format": dict(
            help=HELP_INFO["gene_format"], default='txt',
            choices=['txt', 'bed6', 'bed12'],
        ),
        "exclude_supplementary": dict(
            help=HELP_INFO["exclude_supplementary"], action="store_true"
        ),
        "exclude_duplicate": dict(
            help=HELP_INFO["exclude_duplicate"], action="store_true"
        ),
        "exclude_secondary": dict(
            help=HELP_INFO["exclude_secondary"], action="store_true"
        ),
        "proper_pair_only": dict(
            help=HELP_INFO["proper_pair_only"], action="store_true"
        ),
    }
    parser.add_argument("bamFile", **arguments["bam"])
    parser.add_argument("-v",   "--version",       **arguments["v"])
    parser.add_argument("-bed", "--bedFile",       **arguments["bed"])
    parser.add_argument("-g",   "--genetable",     **arguments["g"])
    parser.add_argument("-r",   "--chromosome",    **arguments["r"])
    parser.add_argument("-O",   "--outPath",       **arguments["O"])
    parser.add_argument("-q",   "--mapq",          **arguments["q"])
    parser.add_argument("--window_size",           **arguments["window_size"])
    parser.add_argument("--gene_format",           **arguments["gene_format"])
    parser.add_argument("--exclude_supplementary", **arguments["exclude_supplementary"])
    parser.add_argument("--exclude_duplicate",     **arguments["exclude_duplicate"])
    parser.add_argument("--exclude_secondary",     **arguments["exclude_secondary"])
    parser.add_argument("--proper_pair_only",      **arguments["proper_pair_only"])
    return parser

if __name__ == "__main__":
    main()
