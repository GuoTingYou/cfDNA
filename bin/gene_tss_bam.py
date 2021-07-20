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
    BAM, VCF,
    GeneTable,
    PassSnpSelector,
    ReadSNP,
    # function
    argParse_template,
    timemain,
    outfile,
)
HELP_INFO = {
    "bam": "Input: <cfDNA.bam> with absolute path.",
    "bg":  "Input: <background.VCF> with absolute path. e.g. maternal",
    "it":  "Input: <interested.VCF> with absolute path. e.g. fetal",
    "g":   "Input: <tissue_specific_gene> with absolute path.",
    "r":   "Chromosome to compare VCF.",
    "O":   "Output: Absolute path for output file.",
    "q":   "Parameter: Mapping Quality of READ in BAM above `q` will be pileuped.",
    "Q":   "Parameter: Base Quality of READ in BAM above `Q` will be pileuped.",
    "extend": "Parameter:",
    "gene_format": "Parameter:",
    "exclude_supplementary": "Parameter:",
    "exclude_duplicate": "Parameter:",
    "exclude_secondary": "Parameter:",
    "proper_pair_only": "Parameter:",
}
@timemain
def main():

    ### I/O ###

    ioargs = argParse(4)

    # input
    CHROM = ioargs.chromosome
    bamFile = BAM(ioargs.bamFile)
    genetable = GeneTable(ioargs.genetable, fmt=ioargs.gene_format)

    # output
    outBamBG = outfile(
        ioargs.bamFile, outPath = ioargs.outPath,
        suffix = f'cfDNA.bg_specific.{CHROM}.{tissue(ioargs.genetable)}.tss.bam'
    )
    outBamBS = outfile(
        ioargs.bamFile, outPath = ioargs.outPath,
        suffix = f'cfDNA.bg_share.{CHROM}.{tissue(ioargs.genetable)}.tss.bam'
    )
    outBamIT = outfile(
        ioargs.bamFile, outPath = ioargs.outPath,
        suffix = f'cfDNA.it_specific.{CHROM}.{tissue(ioargs.genetable)}.tss.bam'
    )
    outBamIS = outfile(
        ioargs.bamFile, outPath = ioargs.outPath,
        suffix = f'cfDNA.it_share.{CHROM}.{tissue(ioargs.genetable)}.tss.bam'
    )
    outTssBG = outfile(
        ioargs.bamFile, outPath = ioargs.outPath,
        suffix = f'bg_specific.{tissue(ioargs.genetable)}.snp.gz'
    )
    outTssIT = outfile(
        ioargs.bamFile, outPath = ioargs.outPath,
        suffix = f'it_specific.{tissue(ioargs.genetable)}.snp.gz'
    )
    global HEADER
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

    chromosomes = tuple(f'chr{N}' for N in [*range(1, 23), 'X'])

    tss_regions = defaultdict(list)
    for gene in genetable:
        chrom, start, end = gene.tss_region(extend=ioargs.extend)
        tss_regions[chrom.symbol].append(start)
        tss_regions[chrom.symbol].append(end)

    aaab_snps, abaa_snps = [], []

    with pysam.AlignmentFile(outBamBG, 'wb', template=bamFile.bamFile) as bamBG, \
         pysam.AlignmentFile(outBamBS, 'wb', template=bamFile.bamFile) as bamBS, \
         pysam.AlignmentFile(outBamIT, 'wb', template=bamFile.bamFile) as bamIT, \
         pysam.AlignmentFile(outBamIS, 'wb', template=bamFile.bamFile) as bamIS, \
         gzip.open(outTssBG, 'wt') as tssBG, gzip.open(outTssIT, 'wt') as tssIT:

        if CHROM == 'genome':
            for chrom in chromosomes:
                regions = tss_regions[chrom]
                vcfBG = VCF(ioargs.backgroundVCF).fetch(chrom)
                vcfIT = VCF(ioargs.interestedVCF).fetch(chrom)
                find_snps_in_tss(vcfBG, vcfIT, regions, aaab_snps, abaa_snps)
        else:
            regions = tss_regions[CHROM]
            vcfBG = VCF(ioargs.backgroundVCF).fetch(CHROM)
            vcfIT = VCF(ioargs.interestedVCF).fetch(CHROM)
            find_snps_in_tss(vcfBG, vcfIT, regions, aaab_snps, abaa_snps)

        abaa_snps = ReadSNP(abaa_snps)
        aaab_snps = ReadSNP(aaab_snps)
        extract_alignments(bamFile, abaa_snps, (bamBS, bamBG), tssBG, **filters)
        extract_alignments(bamFile, aaab_snps, (bamIS, bamIT), tssIT, **filters)
    sort_index(outBamBG)
    sort_index(outBamBS)
    sort_index(outBamIT)
    sort_index(outBamIS)

def find_snps_in_tss(vcfBG, vcfIT, regions, aaab_snps, abaa_snps):
    for snp in PassSnpSelector(vcfBG, vcfIT):
        if is_intersect(regions, snp):
            if snp.genotype.startswith('AA-ab'):
                aaab_snps.append(snp)
            elif snp.genotype.startswith('AB-aa'):
                abaa_snps.append(snp)

def extract_alignments(bamFile, snpFile, outBams, outTss, **filters):
    shareBam, specialBam = outBams
    print('#{}'.format('\t'.join(HEADER)), file=outTss)

    for pileup, snp in bamFile.pileup(locFile=snpFile, **filters):
        fwdDepthA, rvsDepthA = pileup.depth_of(snp.allele_A) # fwd: forward strand
        fwdDepthB, rvsDepthB = pileup.depth_of(snp.allele_B) # rvs: reverse strand

        insertsize_A = ','.join(pileup.insertsizes_of(snp.allele_A).astype(str))
        insertsize_B = ','.join(pileup.insertsizes_of(snp.allele_B).astype(str))

        DepthA, DepthB = rvsDepthA + fwdDepthA, rvsDepthB + fwdDepthB

        snp = snp.apply(fwdDepthA = fwdDepthA, rvsDepthA = rvsDepthA,
                        fwdDepthB = fwdDepthB, rvsDepthB = rvsDepthB,
                        fraction = snp.estimate_fraction(DepthA, DepthB))

        snp.print(insertsize_A, insertsize_B, fields=HEADER[:11],
                  sep='\t', file=outTss)

        alleles = (snp.allele_A, snp.allele_B)
        share, special = pileup.get_share_special_alignments(alleles)
        for read in share:
            shareBam.write(read)
        for read in special:
            specialBam.write(read)

def sort_index(bamFile: str):
    d, f = os.path.split(bamFile)
    tmpFile = os.path.join(d, f'tmp.{f}')
    pysam.sort("-o", tmpFile, bamFile)
    pysam.index(tmpFile)
    os.system(f"mv {tmpFile} {bamFile}")
    os.system(f"mv {tmpFile}.bai {bamFile}.bai")

def is_intersect(tss_regions, snp):
    return bisect.bisect(tss_regions, snp.start) % 2 == 1

def tissue(genetable: str):
    return os.path.basename(genetable).split('.')[0]

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
        "bam": dict(help=HELP_INFO["bam"], required=True),
        "bg":  dict(help=HELP_INFO["bg"],  required=True),
        "it":  dict(help=HELP_INFO["it"],  required=True),
        "g":   dict(help=HELP_INFO["g"],   required=True),
        "r":   dict(help=HELP_INFO["r"], default='genome'),
        "O":   dict(help=HELP_INFO["O"], default=os.getcwd()),
        "q":   dict(help=HELP_INFO["q"], type=int, default=30),
        "Q":   dict(help=HELP_INFO["Q"], type=int, default=20),
        "extend": dict(help=HELP_INFO["extend"], type=int, default=1000),
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
    parser.add_argument("-v",   "--version",       **arguments["v"])
    parser.add_argument("-bam", "--bamFile",       **arguments["bam"])
    parser.add_argument("-bg",  "--backgroundVCF", **arguments["bg"])
    parser.add_argument("-it",  "--interestedVCF", **arguments["it"])
    parser.add_argument("-g",   "--genetable",     **arguments["g"])
    parser.add_argument("-r",   "--chromosome",    **arguments["r"])
    parser.add_argument("-O",   "--outPath",       **arguments["O"])
    parser.add_argument("-q",   "--mapq",          **arguments["q"])
    parser.add_argument("-Q",   "--baq",           **arguments["Q"])
    parser.add_argument("--extend",                **arguments["extend"])
    parser.add_argument("--gene_format",           **arguments["gene_format"])
    parser.add_argument("--exclude_supplementary", **arguments["exclude_supplementary"])
    parser.add_argument("--exclude_duplicate",     **arguments["exclude_duplicate"])
    parser.add_argument("--exclude_secondary",     **arguments["exclude_secondary"])
    parser.add_argument("--proper_pair_only",      **arguments["proper_pair_only"])
    return parser

if __name__ == "__main__":
    main()
