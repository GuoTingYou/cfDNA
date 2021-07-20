import os
import sys
import time
import argparse
from collections import defaultdict
from numpy import array

from snp_based_method import ReadFraction

def argParse(N):
    usage = (f"Usage: python {sys.argv[0]}"
              "-vcfM maternal.vcf,gz -vcfF fetal.vcf.gz -r chrN [opetions...]")
    description = (f"Filename: {sys.argv[0]}\n"
                   f"-*- Author: GuoTingYou -*-\n"
                   f"-*- Version: youyou.tortoise 0.0.0 -*-\n"
                   f"-*- Coding: utf-8 -*-\n")
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    help_info = {
        "abaa": "Input: <fetal_fraction.chrN.abaa.gz> file.",
        "aaab": "Input: <fetal_fraction.chrN.aaab.gz> file.",
        "O": "Output: Absolute output path.",
        "o": "Output: Prefix of output name.",
        "d": "Parameter: Depth of the SNP above `d` will be considered.",
        "w": "Parameter: Window size.",
    }
    arguments = {
        "v": dict(action = "version", version = "youyou.tortoise 0.0.0"),
        "abaa": dict(help=help_info["abaa"], required=True),
        "aaab": dict(help=help_info["aaab"], required=True),
        "O": dict(help=help_info["O"], default=os.getcwd()),
        "o": dict(help=help_info["o"], default='output'),
        "d": dict(help=help_info["d"], type=int, default=10),
        "w": dict(help=help_info["w"], type=int, default=1e5),
    }
    parser.add_argument("-v",    "--version",     **arguments["v"])
    parser.add_argument("-abaa", "--input_ABaa",  **arguments["abaa"])
    parser.add_argument("-aaab", "--input_AAab",  **arguments["aaab"])
    parser.add_argument("-O",    "--out_path",    **arguments["O"])
    parser.add_argument("-o",    "--out_name",    **arguments["o"])
    parser.add_argument("-d",    "--depth",       **arguments["d"])
    parser.add_argument("-w",    "--window_size", **arguments["w"])
    if len(sys.argv) <= N: parser.print_help(); raise IndexError(usage)
    return parser.parse_args()

### Class and Function ###

DEPTH = 10
def _read(file, ReadFraction=ReadFraction):
    global DEPTH
    return (line
        for line in ReadFraction(file)
        if line.check_fraction() and line.check_depth(at_least=DEPTH)
    )

def _next(file):
    try:
        return next(file)
    except StopIteration:
        return ReadFraction.stop()

def sorted_merged_fraction_file(fractionFileABaa, fractionFileAAab):
    fractionFileABaa, fractionFileAAab = _read(fractionFileABaa), _read(fractionFileAAab)
    abaa_site, aaab_site = _next(fractionFileABaa), _next(fractionFileAAab)
    while abaa_site.chrom or aaab_site.chrom:

        while aaab_site.chrom is not None and aaab_site < abaa_site:
            yield aaab_site
            aaab_site = _next(fractionFileAAab)

        while abaa_site.chrom is not None and abaa_site < aaab_site:
            yield abaa_site
            abaa_site = _next(fractionFileABaa)

        while abaa_site.chrom is None and aaab_site.chrom is not None:
            yield aaab_site
            aaab_site = _next(fractionFileAAab)

        while aaab_site.chrom is None and abaa_site.chrom is not None:
            yield abaa_site
            abaa_site = _next(fractionFileABaa)

### Main ###

if __name__ == "__main__":

    print(os.path.realpath(__file__))
    print("\n{:-^81s}".format(
          " Task Start At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))

    ioargs = argParse(2)
    # parameter
    DEPTH = ioargs.depth
    WINDOW = ioargs.window_size
    # required
    fractionFileAAab = ioargs.input_AAab
    fractionFileABaa = ioargs.input_ABaa
    # output
    regional_fraction = defaultdict(list)
    regional_depth = defaultdict(list)

    for site in sorted_merged_fraction_file(fractionFileABaa, fractionFileAAab):
        regionID = site.get_regionID(WINDOW)
        regional_fraction[(regionID, regionID+1)].append(site.fraction)
        regional_depth[(regionID, regionID+1)].append(site.depth_sum())

    with open(ioargs.out_name, 'wt') as f:
        for (wstart, wend), fractions in regional_fraction.items():
            depths = array(regional_depth[(wstart, wend)])
            fractions = array(fractions)
            weighted_mean = (fractions * depths).sum() / depths.sum()
            print(wstart, wend, len(fractions), weighted_mean, sep='\t', file=f)

    print("{:-^81s}\n".format(
          " Task End At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))
    print("Practice_makes_perfect")
