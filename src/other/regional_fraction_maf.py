import os
import sys
import time
import argparse
import gzip
import numpy as np

from snp_based_method import ReadSNP
from maf_based_method import ReadMAF, MafDistribution
from partitioner import RegionPartitioner

def argParse(N):
    usage = (f"Usage: python {sys.argv[0]}"
              "-i <cfDNA.chrN.maf.gz> [opetions...]")
    description = (f"Filename: {sys.argv[0]}\n"
                   f"-*- Author: GuoTingYou -*-\n"
                   f"-*- Version: youyou.tortoise 0.0.0 -*-\n"
                   f"-*- Coding: utf-8 -*-\n")
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    help_info = {
        "maf":"Input: MAF estimation regional fetal fraction.",
        "abaa": "Input: <fetal_fraction.chrN.abaa.gz> file.",
        "aaab": "Input: <fetal_fraction.chrN.aaab.gz> file.",
        "O":  "Absolute output path.",
        "c":  "Criterion for RegionPartitioner.",
        "n":  "Number of pileup sites in a window above n will be consider.",
        "bw": "Parameter pass to scipy.stats.gaussian_kde bw_method.",
        "compare": "If given, compare with vcf_based_method regional fraction.",
        "skip_indel": "If given, skip indel while estimate maf distribution.",
    }
    arguments = {
        "v": dict(action = "version", version = "youyou.tortoise 0.0.0"),
        "maf": dict(help=help_info["maf"],  required=True),
        "abaa":dict(help=help_info["abaa"], default=None),
        "aaab":dict(help=help_info["aaab"], default=None),
        "O": dict(help=help_info["O"], default=os.getcwd()),
        "c": dict(help=help_info["c"], default='gene:1000'),
        "n": dict(help=help_info["n"], type=int, default=0),
        "bw":dict(help=help_info["bw"], type=float, default=.05),
        "compare":   dict(help=help_info["compare"], action="store_true"),
        "skip_indel":dict(help=help_info["skip_indel"], action="store_true"),
    }
    parser.add_argument("-v", "--version",      **arguments["v"])
    parser.add_argument("-maf", "--input_maf",  **arguments["maf"])
    parser.add_argument("-abaa","--input_ABaa", **arguments["abaa"])
    parser.add_argument("-aaab","--input_AAab", **arguments["aaab"])
    parser.add_argument("-O", "--out_path",     **arguments["O"])
    parser.add_argument("-c", "--criterion",    **arguments["c"])
    parser.add_argument("-n", "--num_pups",     **arguments["n"])
    parser.add_argument("-bw", "--bin_width",   **arguments["bw"])
    parser.add_argument("--compare",            **arguments["compare"])
    parser.add_argument("--skip_indel",         **arguments["skip_indel"])
    if len(sys.argv) <= N:
        parser.print_help()
        raise IndexError("Not enough arguments were given.")
    return parser.parse_args()

def outfile(o_path, o_name, region_size):
    o_path = o_path.rstrip('/')
    o_name = os.path.basename(o_name).replace('maf.gz', f'rff.{region_size}.gz')
    return os.path.join(o_path, o_name)

### Main ###

if __name__ == "__main__":

    print(os.path.realpath(__file__))
    print("\n{:-^81s}".format(
          " Task Start At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))

    ### I/O ###

    ioargs = argParse(2)

    # parameter

    SKIP_INDEL = ioargs.skip_indel
    NPUPS      = ioargs.num_pups
    BIN_WIDTH  = ioargs.bin_width

    if ioargs.compare:
        fractionFileAAab = ioargs.input_AAab
        fractionFileABaa = ioargs.input_ABaa
        if fractionFileAAab is None or fractionFileABaa is None:
            raise SyntaxError("If `--compare` is given, both "
                              "`-abaa` and `-aaab` should be given.")
        else:
            from regional_fraction_snp import sorted_merged_fraction_file
            vcf_based_fraction_file = (i
                for i in sorted_merged_fraction_file(fractionFileABaa, fractionFileAAab)
            )

    # required

    mafFile = ReadMAF(ioargs.input_maf)
    criterion, requirement = ioargs.criterion.split(':', 1)
    partition = RegionPartitioner(criterion=criterion, requirement=requirement)

    # output

    o_file = outfile(ioargs.out_path, ioargs.input_maf, criterion)

    with gzip.open(o_file, 'wt') as f:
        print(('#chrom start end region_size depth nmafs '
               'AAab_fraction ABaa_fraction fraction').replace(' ', '\t'), file=f)

        graph = True

        for region in partition(mafFile, skip_indel=SKIP_INDEL):
            if len(region.values) < NPUPS: continue

            maf_distribution = MafDistribution(
                region.depths, region.values, bw_method=BIN_WIDTH)

            maf_distribution.find_peaks(prominence=0.1) # After `find_peaks` method is called,
                                                        # self.peaks and self.distribution
                                                        # attributes are setted.
            if graph:
                outfig = f'maf_distribution.{region.chrom}-{region.start}-{region.end}.pdf'
                maf_distribution.graph(outfig=outfig) # `graph` method should be called after
                                                      # `self.distribution` has been setted.
                graph = False

            # `estimate_fraction` method should be called after `self.peaks` has been setted.
            aaab_fraction, abaa_fraction, weighted_fraction = maf_distribution.estimate_fraction()

            print(region.chrom, region.start, region.end, region.size,
                  sum(region.depths), maf_distribution.nloci,
                  aaab_fraction, abaa_fraction, weighted_fraction, sep='\t', file=f)

    if ioargs.compare:

        with gzip.open(f'{ioargs.out_path}/regional_fraction.{chrom}.{criterion}-{requirement}.gz', 'wt') as f:
            print('\t'.join('#chrom start end region_size depth nloci '
                            'fraction'.split()), file=f)

            for region in partition(vcf_based_fraction_file, skip_indel=SKIP_INDEL):

                depths = np.array(region.depths)
                values = np.array(region.values)
                weighted_fraction = values.dot(depths / depths.sum()).mean()

                print(region.chrom, region.start, region.end, region.size,
                      sum(region.depths), len(region.values), weighted_fraction, sep='\t', file=f)

    print("{:-^81s}\n".format(
          " Task End At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))
    print("Practice_makes_perfect")
