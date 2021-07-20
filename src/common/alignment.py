import pysam

from ..base import CellFreeDNA
from ..base import PysamPileupRead
from ..base import PysamPileupColumn
from ..util import pipeline

from .location import FetchLocation, PileupLocation
from .pileup import Pileup

class BAM:

    def __init__(self, filename):
        self.bamFile = pysam.AlignmentFile(filename)

    def __iter__(self):
        return map(CellFreeDNA, self.bamFile)

    def fetch(self, *args, locFile=None, is_region=False, **kwargs):
        if locFile:
            for params, loc in FetchLocation(locFile):
                params.update(kwargs)
                if is_region:
                    yield self._fetch(**params), loc
                else:
                    for read in self._fetch(**params):
                        yield read, loc
        else:
            yield from self._fetch(*args, **kwargs)

    def _fetch(self, *args, MQ=0,
               exclude_supplementary=True,
               exclude_duplicate=True,
               exclude_secondary=True,
               proper_pair_only=True,
               unique_score=None, mismatch=None, **kwargs):

        reads_options = {
            'MQ': MQ,
            'exclude_supplementary': exclude_supplementary,
            'exclude_duplicate': exclude_duplicate,
            'exclude_secondary': exclude_secondary,
            'proper_pair_only': proper_pair_only,
            'unique_score': unique_score,
            'mismatch': mismatch,
        }
        for read in map(CellFreeDNA, self.bamFile.fetch(*args, **kwargs)):
            if read.apply_filters(**reads_options):
                yield read

    def pileup(self, *args, locFile=None, is_region=False, **kwargs):
        if locFile:
            for params, loc in PileupLocation(locFile):
                params.update(kwargs)
                if is_region:
                    yield self._pileup(**params), loc
                else:
                    for pileup in self._pileup(**params):
                        yield pileup, loc
        else:
            yield from self._pileup(*args, **kwargs)

    def _pileup(self, *args, BQ=13, MQ=0,
                exclude_ins=True, exclude_del=True,
                exclude_supplementary=True,
                exclude_duplicate=True,
                exclude_secondary=True,
                proper_pair_only=True,
                unique_score=None, mismatch=None, **kwargs):

        indel_options = {
            'exclude_ins': exclude_ins,
            'exclude_del': exclude_del,
        }
        reads_options = {
            'MQ': MQ,
            'exclude_supplementary': exclude_supplementary,
            'exclude_duplicate': exclude_duplicate,
            'exclude_secondary': exclude_secondary,
            'proper_pair_only': proper_pair_only,
            'unique_score': unique_score,
            'mismatch': mismatch,
        }
        for pileupcolumn in map(PysamPileupColumn,
                                self.bamFile.pileup(*args, **kwargs)):
            pileupcolumn.apply_filters(BQ=BQ)
            pileupreads = (
                PysamPileupRead(p).apply_filters(**indel_options)
                for p in pileupcolumn.pileups
            )
            reads = (
                CellFreeDNA(p.alignment).apply_filters(**reads_options)
                for p in pileupcolumn.pileups
            )
            inclusive_reads = [
               pileupread_is_pass and read_is_pass
               for pileupread_is_pass, read_is_pass in zip(pileupreads, reads)
            ]
            yield Pileup(pileupcolumn, inclusive_reads)
