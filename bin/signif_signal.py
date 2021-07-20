import os
import sys
import gzip
import argparse
import numpy as np
import pyBigWig as pbw

from os.path import realpath, abspath, dirname, join
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

from src import ReadRegionalFeature, GeneTable, ReadBed
from src import argParse_template, outfile

USAGE = f"python {sys.argv[0]}"

HELP_INFO = {
    "i": "Input",
    "p": "Input",
    "w": "Input",
    "O": "Output",
    "o": "Output",
    "fdr": "Parameter",
}
def main():
    global SIGNAL_PLACENTA, SIGNAL_PBT_CELL, THRESHOLD_FDR

    ioargs = argParse(2)
    rfzFile = ioargs.regionFile
    SIGNAL_PBT_CELL = ioargs.wbc       # bigwig file
    SIGNAL_PLACENTA = ioargs.placenta  # bigwig file
    THRESHOLD_FDR = ioargs.threshold_fdr
    outFile = outfile(rfzFile, outPath=ioargs.outPath,
                      prefix=ioargs.outName, suffix='acc.gz')

    with gzip.open(outFile, 'wt') as out:
        #for region in ReadRegionalFeature(rfzFile):
        #for region in GeneTable(rfzFile):
        for region in ReadBed(rfzFile):
            arr_placenta = region.signal_from_bigwig(SIGNAL_PLACENTA)
            arr_pbt_cell = region.signal_from_bigwig(SIGNAL_PBT_CELL)
            arr_placenta = arr_placenta[~np.isnan(arr_placenta)]
            arr_pbt_cell = arr_pbt_cell[~np.isnan(arr_pbt_cell)]
            if arr_placenta.size == 0 or arr_pbt_cell.size == 0:
                continue

            # accessible region
            acc_region_placenta = number_of_signif_locus(arr_placenta)
            acc_region_pbt_cell = number_of_signif_locus(arr_pbt_cell)

            # accessibility score
            acc_scores_placenta = sum_of_signif_p_values(arr_placenta)
            acc_scores_pbt_cell = sum_of_signif_p_values(arr_pbt_cell)

            # accessible fraction
            acc_fraction_placenta = acc_region_placenta / region.size
            acc_fraction_pbt_cell = acc_region_pbt_cell / region.size

            # output
            #region.print(*region.tss_region(),
            region.print(acc_region_placenta,   acc_region_pbt_cell,
                         acc_scores_placenta,   acc_scores_pbt_cell,
                         acc_fraction_placenta, acc_fraction_pbt_cell,
                         sep='\t', file=out)

def number_of_signif_locus(arr):
    global THRESHOLD_FDR
    threshold = -np.log10(THRESHOLD_FDR / arr.size)
    return np.where(arr >= threshold, 1, 0).sum()

def sum_of_signif_p_values(arr):
    global THRESHOLD_FDR
    threshold = -np.log10(THRESHOLD_FDR / arr.size)
    return np.where(arr >= threshold, arr, 0).sum()

@argParse_template
def argParse(N):
    global USAGE, HELP_INFO
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = USAGE, description = __doc__,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    arguments = {
        "v": dict(action="version", version="youyou 0.0.0"),
        "i": dict(help=HELP_INFO["i"]),
        "p": dict(help=HELP_INFO["p"], required=True, type=pbw.open),
        "w": dict(help=HELP_INFO["w"], required=True, type=pbw.open),
        "O": dict(help=HELP_INFO["O"], default=os.getcwd()),
        "o": dict(help=HELP_INFO["o"], default=1),
        "fdr": dict(help=HELP_INFO["fdr"], type=float, default=0.01),
    }
    parser.add_argument("regionFile",   **arguments["i"])
    parser.add_argument("--placenta",   **arguments["p"])
    parser.add_argument("--wbc",        **arguments["w"])
    parser.add_argument("-v", "--version", **arguments["v"])
    parser.add_argument("-O", "--outPath", **arguments["O"])
    parser.add_argument("-o", "--outName", **arguments["o"])
    parser.add_argument("-fdr", "--threshold_fdr", **arguments["fdr"])
    return parser

if __name__ == "__main__":
    main()
