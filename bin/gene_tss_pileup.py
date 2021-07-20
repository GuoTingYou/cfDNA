from os.path import realpath, abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

import argparse
import bisect
import gzip
import pysam
from collections import defaultdict, namedtuple

from src import (
    # class
    GeneTable,
    # function
    argParse_template,
    timemain,
    outfile,
)
HELP_INFO = {
    "p": "Input: <pileup.gz>.",
    "g": "Input: <tissue_specific_gene> with absolute path.",
    "r": "Chromosome to compare VCF.",
    "O": "Output: Absolute path for output file.",
    "q": "Parameter: Mapping Quality of READ in BAM above `q` will be pileuped.",
    "Q": "Parameter: Base Quality of READ in BAM above `Q` will be pileuped.",
    "extend": "Parameter:",
    "gene_format": "Parameter:",
}
@timemain
def main():

    ### I/O ###

    ioargs = argParse(1)

    # input
    CHROM = ioargs.chromosome
    pupFile = ioargs.pileupFile
    genetable = GeneTable(ioargs.genetable, fmt=ioargs.gene_format)

    # output
    outFile = outfile(
        ioargs.pileupFile, outPath = ioargs.outPath,
        prefix = 3, suffix = 'tss_coverage.bed.gz',
    )
    global HEADER
    HEADER = [
        'chrom',  # Chromosome of the SNP (with chr prefix).
        'start',  # Start position relative to TSS (0-based coordinate).
        'end',    # End position relative to TSS (end exclusive).
        'symbol', # Symbol of Gene.
        'depth',  # Depth of coverage.
        'tss',    # Transcription start site of gene.
    ]

    ### main ###

    chromosomes = tuple(f'chr{N}' for N in [*range(1, 23), 'X'])

    global TSS_REGIONS, TSS_GENES
    TSS_REGIONS = defaultdict(list)
    TSS_GENES = defaultdict(list)

    for gene in genetable:
        chrom, start, end = gene.tss_region(extend=ioargs.extend)
        TSS_REGIONS[chrom.symbol].append(start)
        TSS_REGIONS[chrom.symbol].append(end)
        TSS_GENES[chrom.symbol].append(gene)

    with gzip.open(outFile, 'wt') as out:
        print('#{}'.format('\t'.join(HEADER)), file=out)
        for locus in parse_pileup_file(pupFile):
            gene = find_intersected_gene(locus)
            if gene:
                tss = gene.transcription_start_site
                relative_pos = locus.pos - tss
                name = getattr(gene, 'symbol', getattr(gene, 'name'))
                print(gene.chrom.symbol, relative_pos, relative_pos + 1,
                      name, locus.depth, tss, sep='\t', file=out)

def find_intersected_gene(locus):
    index = bisect.bisect(TSS_REGIONS[locus.chrom], locus.pos-1)
    if index % 2 == 1:
        return TSS_GENES[locus.chrom][index // 2]
    else:
        return None

def parse_pileup_file(pupFile):
    # Locus will convert `pos` to 0-based
    Locus = namedtuple('Locus', 'chrom, pos, depth')
    with gzip.open(pupFile, 'rt') as f:
        f = (ln.strip().split('\t') for ln in f)
        for chrom, pos, _, depth, *_ in f:
            yield Locus(chrom, int(pos)-1, int(depth))

@argParse_template
def argParse(N):
    usage = (f"python {sys.argv[0]} "
              "-bam <cfDNA.bam> -bg <background.VCF> -it <interested.VCF> [-O ...]")
    description = f"Filename: {sys.argv[0]}"
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    arguments = {
        "v": dict(action="version", version="youyou 0.0.0"),
        "p": dict(help=HELP_INFO["p"]),
        "g": dict(help=HELP_INFO["g"], required=True),
        "r": dict(help=HELP_INFO["r"], default='genome'),
        "O": dict(help=HELP_INFO["O"], default=os.getcwd()),
        "q": dict(help=HELP_INFO["q"], type=int, default=30),
        "Q": dict(help=HELP_INFO["Q"], type=int, default=20),
        "extend": dict(help=HELP_INFO["extend"], type=int, default=1000),
        "gene_format": dict(
            help=HELP_INFO["gene_format"], default='txt',
            choices=['txt', 'bed6', 'bed12'],
        ),
    }
    parser.add_argument("pileupFile", **arguments["p"])
    parser.add_argument("-v",   "--version",    **arguments["v"])
    parser.add_argument("-g",   "--genetable",  **arguments["g"])
    parser.add_argument("-r",   "--chromosome", **arguments["r"])
    parser.add_argument("-O",   "--outPath",    **arguments["O"])
    parser.add_argument("-q",   "--mapq",       **arguments["q"])
    parser.add_argument("-Q",   "--baq",        **arguments["Q"])
    parser.add_argument("--extend",             **arguments["extend"])
    parser.add_argument("--gene_format",        **arguments["gene_format"])
    return parser

if __name__ == "__main__":
    main()
