from os.path import realpath, abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

import argparse
import gzip
from collections import defaultdict

from src import IntersectFeatureSelector
from src import ReadRegionalFeature
from src import argParse_template
from src import timemain
from src import outfile

USAGE = f"python {sys.argv[0]}"

HELP_INFO = {
    "i": "Input",
    "f": "Input",
    "O": "Output",
    "o": "Output",
    "DNaseHotspot": "Parameter",
    "NarrowPeak":   "Parameter",
    "BroadPeak":    "Parameter",
    "GappedPeak":   "Parameter",
}
@timemain
def main():
    ioargs = argParse(1)
    rfzFile = ReadRegionalFeature(ioargs.input)

    if ioargs.DNaseHotspot:
        from src import ReadDNaseHotspot as ReadFeature
        SUFFIX = 'DNaseHotspot'
    elif ioargs.NarrowPeak:
        from src import ReadNarrowPeak as ReadFeature
        SUFFIX = 'NarrowPeak'
    elif ioargs.BroadPeak:
        from src import ReadBroadPeak as ReadFeature
        SUFFIX = 'BroadPeak'
    elif ioargs.GappedPeak:
        from src import ReadGappedPeak as ReadFeature
        SUFFIX = 'GappedPeak'
    else:
        raise IndexError("Type of feature must be specified.")
    RoadMap = ReadFeature(ioargs.feature)
    RoadMap = ReadFeature((f for f in sorted(RoadMap)))
    outFile = outfile(ioargs.input, outPath = ioargs.outPath,
                      suffix = f'intersect.{SUFFIX}.gz')
    f = open('placenta.DNaseHotspot.highconf.bed', 'w')
    for r in RoadMap:
        r.print(sep='\t', file=f)
    f.close()
    sys.exit()
    mapping = defaultdict(list)

    for region, feature, *fractions in IntersectFeatureSelector(rfzFile, RoadMap):
        if feature is None:
            mapping[region] = []
            continue
        mapping[region].append((feature, *fractions))

    with gzip.open(outFile, 'wt') as out:
        for region, values in mapping.items():
            score = 0
            weights = 0
            overlap = 0

            for feature, _, overlap_fraction in values:
                score += feature.signal_value * overlap_fraction \
                      if hasattr(feature, 'signal_value') \
                    else feature.score * overlap_fraction
                weights += overlap_fraction
                overlap += feature.size * overlap_fraction
            else:
                normalized_feature_count = len(values) / region.size * 1e6
                weighted_signal_value = score / weights
                accessible_fraction = overlap / region.size

            region.print(
                f'{normalized_feature_count:.4f}',
                f'{weighted_signal_value:.4f}',
                f'{accessible_fraction:.4f}',
                sep='\t', file=out)

@argParse_template
def argParse(N):
    global USAGE, HELP_INFO
    parser = argparse.ArgumentParser(
             prog = sys.argv[0], usage = USAGE, description = __doc__,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    arguments = {
        "v": dict(action="version", version="youyou 0.0.0"),
        "i": dict(help=HELP_INFO["i"], default=sys.stdin),
        "f": dict(help=HELP_INFO["f"], required=True),
        "O": dict(help=HELP_INFO["O"], default=os.getcwd()),
        "o": dict(help=HELP_INFO["o"], default='intersect'),
        "DNaseHotspot": dict(help=HELP_INFO["DNaseHotspot"], action="store_true"),
        "NarrowPeak":   dict(help=HELP_INFO["NarrowPeak"], action="store_true"),
        "BroadPeak":    dict(help=HELP_INFO["BroadPeak"], action="store_true"),
        "GappedPeak":   dict(help=HELP_INFO["GappedPeak"], action="store_true"),
    }
    parser.add_argument("-v", "--version", **arguments["v"])
    parser.add_argument("-i", "--input",   **arguments["i"])
    parser.add_argument("-f", "--feature", **arguments["f"])
    parser.add_argument("-O", "--outPath", **arguments["O"])
    parser.add_argument("-o", "--outName", **arguments["o"])
    parser.add_argument("--DNaseHotspot",  **arguments["DNaseHotspot"])
    parser.add_argument("--NarrowPeak",    **arguments["NarrowPeak"])
    parser.add_argument("--BroadPeak",     **arguments["BroadPeak"])
    parser.add_argument("--GappedPeak",    **arguments["GappedPeak"])
    return parser

if __name__ == "__main__":
    main()
