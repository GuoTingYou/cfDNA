import os
import sys
import gzip
import argparse
import numpy as np
import pandas as pd

from os.path import realpath, abspath, dirname, join
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

from src import ReadRegionalFeature
from src import argParse_template, outfile

USAGE = f"python {sys.argv[0]}"

HELP_INFO = {
    "i": "Input",
    "O": "Output",
    "o": "Output",
    "fdr": "Parameter",
}
def main():
    global SIGNAL_PLACENTA, SIGNAL_PBT_CELL, THRESHOLD_FDR

    ioargs = argParse(0)
    rfzFile = ioargs.regionFile
    THRESHOLD_FDR = ioargs.threshold_fdr
    outFile = outfile(rfzFile, outPath=ioargs.outPath,
                      prefix=4, suffix='z.gz')

    index_constant, fractions_constant = [], []
    idnex_fetal, fractions_fetal = [], []
    index_maternal, fractions_maternal = [], []

    with gzip.open(outFile, 'wt') as out:
        for region in ReadRegionalFeature(rfzFile):
            if region.name.startswith('constant'):
                index_constant.append(region)
                fractions_constant.append(abs(region.fraction))
            elif region.name.startswith('fetal'):
                idnex_fetal.append(region)
                fractions_fetal.append(abs(region.fraction))
            elif region.name.startswith('maternal'):
                index_maternal.append(region)
                fractions_maternal.append(abs(region.fraction))
        else:
            fractions_constant = pd.Series(fractions_constant, index=index_constant)
            fractions_fetal = pd.Series(fractions_fetal, index=idnex_fetal)
            fractions_maternal = pd.Series(fractions_maternal, index=index_maternal)

            # drop extreme outlier
            median, std = fractions_constant.quantile(0.5), fractions_constant.std()
            outliers = (
                (fractions_constant > median + 1.5 * std) |
                (fractions_constant < median - 1.5 * std)
            )
            fractions_constant = fractions_constant[~outliers]

            # z scores
            median, std = fractions_constant.quantile(0.5), fractions_constant.std()
            fractions_fetal = (fractions_fetal - median) / std
            fractions_maternal = (fractions_maternal - median) / std

        regions = pd.concat([
            fractions_constant,
            fractions_fetal,
            fractions_maternal,
        ]).sort_index()
        for region in regions.index:
            zscore = regions[region]
            print(*region, zscore, sep='\t', file=out)

@argParse_template
def argParse(N):
    global USAGE, HELP_INFO
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = USAGE, description = __doc__,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    arguments = {
        "v": dict(action="version", version="youyou 0.0.0"),
        "i": dict(help=HELP_INFO["i"]),
        "O": dict(help=HELP_INFO["O"], default=os.getcwd()),
        "o": dict(help=HELP_INFO["o"], default=1),
        "fdr": dict(help=HELP_INFO["fdr"], type=float, default=0.01),
    }
    parser.add_argument("regionFile",   **arguments["i"])
    parser.add_argument("-v", "--version", **arguments["v"])
    parser.add_argument("-O", "--outPath", **arguments["O"])
    parser.add_argument("-o", "--outName", **arguments["o"])
    parser.add_argument("-fdr", "--threshold_fdr", **arguments["fdr"])
    return parser

if __name__ == "__main__":
    main()
