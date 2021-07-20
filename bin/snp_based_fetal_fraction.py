from os.path import realpath, abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

import argparse
import gzip

from src import (
    # class
    BAM, VCF,
    PassSnpSelector,
    # function
    argParse_template,
    timemain,
    outfile,
)
@argParse_template
def argParse(N):
    usage = (f"python {sys.argv[0]} "
              "-r chrN -bam <cfDNA.bam> -bg <background.VCF> -it <interested.VCF> [-O ...]")
    description = f"Filename: {sys.argv[0]}"
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    help_info = {
        "r":   "Chromosome to compare VCF.",
        "bam": "Input: <cfDNA.bam> with absolute path.",
        "bg":  "Input: <background.VCF> with absolute path. e.g. maternal",
        "it":  "Input: <interested.VCF> with absolute path. e.g. fetal",
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
        "bg":  dict(help=help_info["bg"],  required=True),
        "it":  dict(help=help_info["it"],  required=True),
        "O":   dict(help=help_info["O"], default=os.getcwd()),
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
    parser.add_argument("-bg",  "--backgroundVCF", **arguments["bg"])
    parser.add_argument("-it",  "--interestedVCF", **arguments["it"])
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

    ioargs = argParse(4)

    # input
    CHROM = ioargs.chromosome
    bamFile = BAM(ioargs.bamFile)
    vcfBG = VCF(ioargs.backgroundVCF).fetch(CHROM)
    vcfIT = VCF(ioargs.interestedVCF).fetch(CHROM)

    # output
    o_fraction = outfile(ioargs.bamFile,
                         outPath=ioargs.outPath,
                         suffix=f'{CHROM}.snp.gz')
    header1 = [
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
    ]
    o_insertsize = outfile(ioargs.bamFile,
                           outPath=ioargs.outPath,
                           suffix=f'insertsize.{CHROM}.gz')
    header2 = [
        'chrom',        # Chromosome of the SNP (with chr prefix).
        'start',        # Start position of the SNP (0-based coordinate).
        'end',          # End position of the SNP (end exclusive).
        'genotype',     # Combination genotype of the SNP locus.
        'allele_A',     # Share allele of combination genotype of the SNP.
        'allele_B',     # Special allele of combination genotype of the SNP.
        'insertsize_A', # Insertsize of cfDNA with share allele.
        'insertsize_B', # Insertsize of cfDNA with special allele.
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

    with gzip.open(o_fraction, 'wt') as f1, gzip.open(o_insertsize, 'wt') as f2:

        # header for out files
        print('#{}'.format('\t'.join(header1)), file=f1)
        print('#{}'.format('\t'.join(header2)), file=f2)

        for pileup, snp in bamFile.pileup(locFile=PassSnpSelector(vcfBG, vcfIT), **filters):

            fwdDepthA, rvsDepthA = pileup.depth_of(snp.allele_A) # fwd: forward strand
            fwdDepthB, rvsDepthB = pileup.depth_of(snp.allele_B) # rvs: reverse strand

            insertsize_A = ','.join(pileup.insertsizes_of(snp.allele_A).astype(str))
            insertsize_B = ','.join(pileup.insertsizes_of(snp.allele_B).astype(str))

            DepthA = rvsDepthA + fwdDepthA
            DepthB = rvsDepthB + fwdDepthB

            snp = snp.apply(fwdDepthA = fwdDepthA, rvsDepthA = rvsDepthA,
                            fwdDepthB = fwdDepthB, rvsDepthB = rvsDepthB,
                            fraction = snp.estimate_fraction(DepthA, DepthB))

            snp.print(fields=header1, sep='\t', file=f1)
            snp.print(insertsize_A, insertsize_B, fields=header2[:5], sep='\t', file=f2)

if __name__ == "__main__":
    main()
