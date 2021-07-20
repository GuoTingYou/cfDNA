from os.path import realpath, abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

import argparse
import gzip

from src import (
    # class
    BAM, VCF,
    PassLocSelector,
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
        "r":   "Chromosome to pileup from location of passed variants in given <VCF>.",
        "bam": "Input: <cfDNA.BAM> with absolute path.",
        "vcf": "Input: <cfDNA.VCF> with absolute path.",
        "O":   "Output: Absolute path for output file.",
        "q":   "Parameter: Mapping Quality of READ in BAM above `q` will be pileuped.",
        "Q":   "Parameter: Base Quality of READ in BAM above `Q` will be pileuped.",
        "exclude_supplementary": "Parameter:",
        "exclude_duplicate": "Parameter:",
        "exclude_secondary": "Parameter:",
        "proper_pair_only": "Parameter:",
    }
    arguments = {
        "v":   dict(action="version", version="youyou 0.0.0"),
        "r":   dict(help=help_info["r"],   required=True),
        "bam": dict(help=help_info["bam"], required=True),
        "vcf": dict(help=help_info["vcf"], required=True),
        "O":   dict(help=help_info["O"],   default=os.getcwd()),
        "q":   dict(help=help_info["q"], type=int, default=30),
        "Q":   dict(help=help_info["Q"], type=int, default=20),
        "exclude_supplementary": dict(
            help=help_info["exclude_supplementary"], action="store_true"),
        "exclude_duplicate": dict(
            help=help_info["exclude_duplicate"], action="store_true"),
        "exclude_secondary": dict(
            help=help_info["exclude_secondary"], action="store_true"),
        "proper_pair_only": dict(
            help=help_info["proper_pair_only"], action="store_true"),
    }
    parser.add_argument("-v",   "--version",       **arguments["v"])
    parser.add_argument("-r",   "--chromosome",    **arguments["r"])
    parser.add_argument("-bam", "--bamFile",       **arguments["bam"])
    parser.add_argument("-vcf", "--vcfFile",       **arguments["vcf"])
    parser.add_argument("-O",   "--outPath",       **arguments["O"])
    parser.add_argument("-q",   "--mapq",          **arguments["q"])
    parser.add_argument("-Q",   "--baq",           **arguments["Q"])
    parser.add_argument("--exclude_supplementary", **arguments["exclude_supplementary"])
    parser.add_argument("--exclude_duplicate",     **arguments["exclude_duplicate"])
    parser.add_argument("--exclude_secondary",     **arguments["exclude_secondary"])
    parser.add_argument("--proper_pair_only",      **arguments["proper_pair_only"])
    return parser

@timemain
def main():

    ### I/O ###

    ioargs = argParse(3)

    # input
    CHROM = ioargs.chromosome
    bamFile = BAM(ioargs.bamFile)
    vcfFile = VCF(ioargs.vcfFile).fetch(CHROM)

    # output
    o_file = outfile(ioargs.bamFile, outPath=ioargs.outPath, suffix=f'{CHROM}.maf.gz')
    header = [
        'chrom',                  # Chromosome of the pileup locus (with chr prefix).
        'start',                  # Start position of the pileup locus (0-based coordinate).
        'end',                    # End position of the pileup locus (end exclusive).
        'depth',                  # Depth of the pileup locus.
        'major_allele',           # The most common allele type of the pileup locus.
        'minor_allele',           # The second common allele type of the pileup locus.
        'fwdDepthA',              # Depth of forward oriented cfDNA with major allele.
        'rvsDepthA',              # Depth of reverse oriented cfDNA with major allele.
        'fwdDepthB',              # Depth of forward oriented cfDNA with minor allele.
        'rvsDepthB',              # Depth of reverse oriented cfDNA with minor allele.
        'major_allele_frequency', # Allele frequency of major allele.
        'minor_allele_frequency', # Allele frequency of minor allele.
    ]
    # parameters
    filters = {
        'BQ': ioargs.baq,  # base alignment quality
        'MQ': ioargs.mapq, # maping quality
        'exclude_supplementary': ioargs.exclude_supplementary, # default False
        'exclude_duplicate': ioargs.exclude_duplicate,         # default False
        'exclude_secondary': ioargs.exclude_secondary,         # default False
        'proper_pair_only': ioargs.proper_pair_only,           # default False
        'unique_score': 10,
        'mismatch': 5,
    }

    ### main ###

    with gzip.open(o_file, 'wt') as f:
        print('#{}'.format("\t".join(header)), file=f)

        for pileup, loc in bamFile.pileup(locFile=PassLocSelector(vcfFile), **filters):

            # After `cal_allele_frequencies` method is called,
            # self.major_allele, self.major_allele_frequency,
            # self.minor_allele, self.minor_allele_frequency,
            # attributes are setted.
            pileup.cal_allele_frequencies(bq_weighted=False)

            fwdDepthA, rvsDepthA = pileup.depth_of(pileup.major_allele) # A for major
            fwdDepthB, rvsDepthB = pileup.depth_of(pileup.minor_allele) # B for minor

            pileup.print(fwdDepthA, rvsDepthA, fwdDepthB, rvsDepthB,
                         pileup.major_allele_frequency,
                         pileup.minor_allele_frequency,
                         fields=header[:6], sep='\t', file=f)

if __name__ == "__main__":
    main()
