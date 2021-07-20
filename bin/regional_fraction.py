from os.path import abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(__file__), '..')))

import argparse
import functools
import gzip

import numpy as np

from src import RegionPartitioner
from src import timemain, argParse_template, outfile

HELP_INFO = {
    "maf": "Input: <chrN.maf.gz> file.",
    "snp": "Input: <chrN.snp.gz> file.",
    "O":  "Output: Absolute path for output file.",
    "o":  "Output: Prefix for output file name.",
    "c":  "Parameter: Criterion for RegionPartitioner.",
    "g":  "Parameter: Gene table that used for defining regions.",
    "n":  "Parameter: Number of loci in a region above `n` will be considered.",
    "bw": "Parameter: Bin width for scipy.stats.gaussian_kde `bw_method` argument.",
    "gene_format": "Parameter: Format of GeneTable.",
    "skip_indel": "Parameter: If given, skip indel while estimate regional fraction.",
    "dependent_partition": (
        "Parameter: If given, use the same partitioner to divide "
        "<mafFile> and <snpFile> into the same regions."
    ),
}
@timemain
def main():

    ### I/O ###

    ioargs = argParse(1)

    if ioargs.mafFile is None and ioargs.snpFile is None:
        raise IndexError('Not enough arguments were given. '
                         'At least `-maf` or `-snp` must be specified.')

    if ioargs.criterion.find(':') == -1:
        raise ValueError(f'Illegal format of "{ioargs.criterion}" for argument '
                          '`-c (--criterion)`, which should be colon (":") separated.')
    else:
        NLOCI      = ioargs.nloci
        BIN_WIDTH  = ioargs.bin_width
        SKIP_INDEL = ioargs.skip_indel
        criterion, requirement = ioargs.criterion.split(':', 1)
        SUFFIX = 'bed' if criterion == 'bed' else ioargs.criterion.replace(":", "_")
        partition = RegionPartitioner(
            criterion = criterion,
            requirement = requirement,
            genetable = ioargs.genetable,
        )
    header = [
        'chrom',    # Chromosome of the defined region (with chr prefix).
        'start',    # Start position of the region (0-based coordinate).
        'end',      # End position of the region (end exclusive).
        'name',     # Name of the defined region (unit: base pair).
        'depth',    # Total depth of the defined region.
        'nloci',    # Number of locus for estimate cfDNA fraction in the region.
        'fraction', # Fraction of cfDNA originated from interested sample.
    ]

    ### Main ###

    if ioargs.mafFile:
        from src import ReadMAF, MafDistribution
        mafFile = ReadMAF(ioargs.mafFile)
        outFile = outfile(ioargs.mafFile,
                          outPath = ioargs.outPath,
                          prefix = ioargs.outName,
                          suffix = f'maf_regional_fraction.{SUFFIX}.gz')

        with gzip.open(outFile, 'wt') as f:
            print('#{}'.format('\t'.join(header)), file=f)

            graph = False

            for region in partition(mafFile, skip_indel=SKIP_INDEL):
                if region.nloci < NLOCI: continue

                maf_distribution = MafDistribution(
                    region.depths, region.values, bw_method=BIN_WIDTH,
                    depth_weighted=True)

                if maf_distribution.kernel is None:
                    continue

                maf_distribution.find_peaks(prominence=0.1)
                # After `find_peaks` method is called,
                # self.peaks and self.distribution
                # attributes are setted.

                if graph:
                    outfig = (f'maf_distribution.{region.chrom}-'
                              f'{region.start}-{region.end}.pdf')
                    # `graph` method should be called after
                    # `self.distribution` has been setted.
                    maf_distribution.graph(outfig=outfig)
                    graph = False

                # `estimate_fraction` method should be called after
                # `self.peaks` has been setted.
                _, _, weighted_fraction = maf_distribution.estimate_fraction()

                region.print(round(maf_distribution.nloci), weighted_fraction,
                             fields=header[:5], sep='\t', file=f)

    if not ioargs.dependent_partition:
        partition.clear()

    if ioargs.snpFile:
        from src import ReadSNP
        snpFile = ReadSNP(ioargs.snpFile)
        outFile = outfile(ioargs.snpFile,
                          outPath = ioargs.outPath,
                          prefix = ioargs.outName,
                          suffix = f'snp_regional_fraction.{SUFFIX}.gz')

        with gzip.open(outFile, 'wt') as f:
            print('#{}'.format('\t'.join(header)), file=f)

            for region in partition(snpFile, skip_indel=SKIP_INDEL):
                if region.nloci < NLOCI: continue

                try:
                    regional_fraction = sum(region.values) / sum(region.depths)
                    # region.values: array of depth of interested allele e.g. fetal
                    # region.depths: array of depth of loci in this region.
                except ZeroDivisionError:
                    continue

                region.print(regional_fraction, fields=header[:6], sep='\t', file=f)

@argParse_template
def argParse(N):
    usage = (f"Usage: python {sys.argv[0]} "
              "-maf <sample.chrN.maf.gz> -snp <sample.chrN.snp.gz> [opetions...]")
    description = f"Filename: {sys.argv[0]}"
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)

    genetable = ('/zfssz2/ST_MCHRI/BIGDATA/USER/guotingyou/'
                 'youyou/Big_Data/ucsc/KnownCanonicalGene.hg19.txt')

    arguments = {
        "v": dict(action = "version", version = "youyou.tortoise 0.0.0"),
        "maf": dict(help=HELP_INFO["maf"], default=None),
        "snp": dict(help=HELP_INFO["snp"], default=None),
        "O": dict(help=HELP_INFO["O"], default=os.getcwd()),
        "o": dict(help=HELP_INFO["o"], default=None),
        "c": dict(help=HELP_INFO["c"], default='dynamic:1000'),
        "g": dict(help=HELP_INFO["g"], default=genetable),
        "n": dict(help=HELP_INFO["n"], type=int, default=0),
        "bw":dict(help=HELP_INFO["bw"], type=float, default=.05),
        "gene_format": dict(
            help=HELP_INFO["gene_format"], default="txt",
            choices=["txt", "bed6", "bed12"],
        ),
        "skip_indel":dict(help=HELP_INFO["skip_indel"], action="store_true"),
        "dependent_partition":dict(help=HELP_INFO["dependent_partition"],
                                   action="store_true"),
    }
    parser.add_argument("-v",   "--version",    **arguments["v"])
    parser.add_argument("-maf", "--mafFile",    **arguments["maf"])
    parser.add_argument("-snp", "--snpFile",    **arguments["snp"])
    parser.add_argument("-O",   "--outPath",    **arguments["O"])
    parser.add_argument("-o",   "--outName",    **arguments["o"])
    parser.add_argument("-c",   "--criterion",  **arguments["c"])
    parser.add_argument("-g",   "--genetable",  **arguments["g"])
    parser.add_argument("-n",   "--nloci",      **arguments["n"])
    parser.add_argument("-bw",  "--bin_width",  **arguments["bw"])
    parser.add_argument("--gene_format",        **arguments["gene_format"])
    parser.add_argument("--skip_indel",         **arguments["skip_indel"])
    parser.add_argument("--dependent_partition",**arguments["dependent_partition"])
    return parser

if __name__ == "__main__":
    main()
