from numpy import r_, zeros, ones
from scipy.signal import fftconvolve

from ..base import GenomeBase
from ..base import GenomeMeta
from ..base import ReadFileBase
from ..base import ReadFileMeta

class ExtraChromosomalCircularDNA(GenomeBase, metaclass = GenomeMeta,
      expand_fields = ('junction_coverage',), type_assert = {'junction_coverage': int}):

    def pileup_param(self):
        flank = int(self.size / 2)
        param = super().pileup_param()
        start = self.start - flank - 1
        start = start if start >= 0 else 0
        end = self.end + flank
        param.update({
            'start': start,
            'stop': end,
        })
        return param

    def estimate_confidence(self, center_depth, flanks_depth):
        #flank = int(self.size / 2)
        #depth_flanks = sum(depth_array * r_[ones(flank), zeros(self.size), ones(flank)])
        #depth_center = sum(depth_array * r_[zeros(flank), ones(self.size), zeros(flank)])
        try:
            return center_depth / flanks_depth
        except ZeroDivisionError:
            return center_depth

    def coverage(self, sequence):
        return 1 - sequence.count('N') / len(sequence)

    def GC_content(self, sequence):
        return (sequence.count('G') + sequence.count('C')) / len(sequence)

    def motif(self, sequence, terminal, N=20):
        if terminal in ("5", "5'"):
            return sequence[:N]
        elif terminal in ("3", "3'"):
            return sequence[-N:]
        else:
            raise ValueError('Valid terminals are ("5", "5\'") or ("3", "3\'")')

class ReadEccDNA(ReadFileBase, metaclass = ReadFileMeta,
                 content_class = ExtraChromosomalCircularDNA):
    def __iter__(self):
        for eccDNA in sorted(super().__iter__()):
            yield eccDNA
