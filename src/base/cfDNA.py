from .genome import GenomeBase, GenomeMeta
from .mixin import AlignedStartMixin

import numpy as np

class CellFreeDNA(AlignedStartMixin):
    """
    General interface class for reads from BAM file.
    """
    def get_cigar_operation_length(self, op):
        return dict(self.cigartuples).get('MIDNSHP=XB'.index(op), 0)

    def get_cigar_operation_ratio(self, op):
        cigar_total_length = self.infer_read_length()
        if cigar_total_length is not None:
            length = self.get_cigar_operation_length(op)
            return length / cigar_total_length
        return None

    def is_uniquely_mapped(self, threshold=None):
        threshold = threshold if threshold else self.query_alignment_length() / 10
        return ((self.best_alignment_score() >= threshold) and
                (self.best_alignment_score() - self.next_alignment_score() > 0))

    def get_edit_distance(self, missing=float('inf')):
        return self.get_tag('NM', missing=missing)

    def best_alignment_score(self, missing=-1):
        return self.get_tag('AS', missing=missing)

    def next_alignment_score(self, missing=-1):
        return self.get_tag('XS', missing=missing)

    def fragment(self, **filter_params):
        if (not self.is_reverse and self.mate_is_reverse and
            not self.is_unmapped and not self.mate_is_unmapped):
            if filter_params:
                if self.apply_filters(**filter_params):
                    return Fragment(self.chrom, self.start, self.end + self.length)
                else:
                    return Fragment.sentinel()
            else:
                return Fragment(self.chrom, self.start, self.end + self.length)
        return Fragment.sentinel()

class Fragment(GenomeBase, metaclass = GenomeMeta):
    def windowed_protection_array(self, window_size: int = 120):
        if window_size % 2 == 1:
            raise ValueError('`window_size` should be an even number.')
        if self.size <= window_size:
            return np.ones(self.size + window_size) * -1
        return np.r_[np.ones(window_size) * -1,
                     np.ones(self.size - window_size),
                     np.ones(window_size) * -1]

    def windowed_protection_coord(self, window_size: int = 120):
        if window_size % 2 == 1:
            raise ValueError('`window_size` should be an even number.')
        return np.arange(self.start - window_size/2,
                         self.end + window_size/2,
                         dtype=np.int32)
