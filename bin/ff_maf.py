from os.path import realpath, abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

import argparse
import gzip

from src import (
    # class
    MafDistribution,
    # function
    argParse_template,
    timemain,
    outfile,
)

@argParse_template
def argParse(N):
    usage = (f"python {sys.argv[0]} "
              "-r chrN -bam <cfDNA.BAM> -vcf <cfDNA.VCF> [-O ...]")
    description = f"Filename: {sys.argv[0]}"
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    help_info = {
        "maf": "Input: <cfDNA.BAM> with absolute path.",
        "l":   "Input: List of mafFile.",
        "O":   "Output: Absolute path for output file.",
        "bw":  "Parameter: Bin Width for non-parameter gaussian kernel."
    }
    arguments = {
        "v":   dict(action="version", version="youyou 0.0.0"),
        "maf": dict(help=help_info["maf"], default=None),
        "l":   dict(help=help_info["l"], type=open, default=None),
        "O":   dict(help=help_info["O"], default=os.getcwd()),
        "bw":  dict(help=help_info["bw"],type=float, default=0.05),
    }
    parser.add_argument("-v",   "--version",   **arguments["v"])
    parser.add_argument("-maf", "--mafFile",   **arguments["maf"])
    parser.add_argument("-l",   "--mafList",   **arguments["l"])
    parser.add_argument("-O",   "--outPath",   **arguments["O"])
    parser.add_argument("-bw",  "--bin_width", **arguments["bw"])
    return parser

@timemain
def main():

    ### I/O ###

    ioargs = argParse(1)
    BIN_WIDTH = ioargs.bin_width

    # input
    mafFile = ioargs.mafFile
    mafList = ioargs.mafList

    if mafFile is None and mafList is None:
        raise IndexError("Either `--mafFile` or `--mafList` must be given.")

    # output
    name = mafList.name if mafList else mafFile
    o_fig = outfile(name, outPath=ioargs.outPath, suffix=f'pdf')

    ### main ###

    depth_list = []
    maf_list = []

    def get_depth_and_maf(line):
        if line.startswith('#'):
            return None, None
        line = line.strip().split('\t')
        depth, maf = line[3], line[-1]
        if depth == "0" or maf == "NA":
            return None, None
        #depth, maf = int(line[2]), float(line[-1])
        return int(depth), float(maf)

    def process_file(mafFile):
        nonlocal depth_list, maf_list
        for depth, maf in map(get_depth_and_maf, mafFile):
            if depth is None:
                continue
            depth_list.append(depth)
            maf_list.append(maf)

    if mafList:
        for mafFile in map(lambda f: gzip.open(f.strip(), 'rt'), mafList):
            process_file(mafFile)

    else:
        process_file(gzip.open(mafFile, 'rt'))

    maf_distribution = MafDistribution(depth_list, maf_list, bw_method=BIN_WIDTH)
    maf_distribution.find_peaks()
    _, _, weighted_fraction = maf_distribution.estimate_fraction()
    print(weighted_fraction)
    maf_distribution.graph(outfig=o_fig)

if __name__ == "__main__":
    main()
