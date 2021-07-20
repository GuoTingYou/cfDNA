class AlignedSegmentFilterMixin:
    """
    A mixin class for PysamAlignedSegmentProxy to apply some basic filterations
    on pysam.AlignedSegment objects.
    """
    def apply_filters(self, *, MQ = 0,
                      exclude_supplementary = False,
                      exclude_duplicate = False,
                      exclude_secondary = False,
                      exclude_qcfail = False,
                      proper_pair_only = False,
                      unique_score = None,
                      mismatch = None):
        return (not self.is_sup(exclude_supplementary) and
                not self.is_dup(exclude_duplicate) and
                not self.is_sec(exclude_secondary) and
                not self.is_qcf(exclude_qcfail) and
                self.is_pp(proper_pair_only) and
                self.is_unique(unique_score) and
                self.is_tolerant_mismatch(mismatch) and
                self.mapping_quality >= MQ)

    def is_sup(self, exclude_supplementary):
        return self.is_supplementary if exclude_supplementary else False

    def is_dup(self, exclude_duplicate):
        return self.is_duplicate if exclude_duplicate else False

    def is_sec(self, exclude_secondary):
        return self.is_secondary if exclude_secondary else False

    def is_qcf(self, exclude_qcfail):
        return self.is_qcfail if exclude_qcfail else False

    def is_pp(self, proper_pair_only):
        return self.is_proper_pair if proper_pair_only else True

    def is_unique(self, unique_score):
        return self.is_uniquely_mapped(unique_score) if unique_score else True

    def is_tolerant_mismatch(self, mismatch):
        return self.get_edit_distance() <= mismatch if mismatch else True

class PileupReadFilterMixin:
    """
    A mixin class for PysamPileupReadProxy to apply some basic filterations
    on pysam.AlignedSegment objects.
    """
    def apply_filters(self, *, exclude_ins=True, exclude_del=True):
        return not self.is_insertion(exclude_ins) and not self.is_deletion(exclude_del)

    def is_insertion(self, exclude_ins):
        return self.is_refskip if exclude_ins else False

    def is_deletion(self, exclude_del):
        return self.is_del if exclude_del else False

class PileupColumnFilterMixin:
    """
    A mixin class for PysamPileupColumnProxy to apply some basic filterations
    on pysam.AlignedSegment objects.
    """
    def apply_filters(self, *, BQ):
        self.set_min_base_quality(BQ)
