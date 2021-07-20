from .filter import (
    AlignedSegmentFilterMixin,
    PileupReadFilterMixin,
    PileupColumnFilterMixin
)

class PysamObjectProxy:
    """
    A proxy class for pysam objects.
    """
    def __init__(self, pysam_object):
        self._pysam_object = pysam_object

    def __str__(self):
        return self._pysam_object.__str__()
    __repr__ = __str__

    def __getattr__(self, attr):
        try:
            return getattr(self._pysam_object, attr)
        except AttributeError:
            raise AttributeError(
                f'Proxy class "{self.__class__.__name__}" or '
                f'Real class "{self._pysam_object.__class__.__name__}" '
                f'object has no attribute or method "{attr}".')

class PysamAlignedSegment(PysamObjectProxy, AlignedSegmentFilterMixin):
    """
    A proxy class for pysam.AlignedSegment.
    """
    def __init__(self, pysam_object):
        super().__init__(pysam_object)
        self.chrom = self.reference_name
        self.start = self.reference_start
        self.end   = self.reference_end
        self.strand = '-' if self.is_reverse else '+'
        self.length = self.template_length
        self.insert_size = self.template_length

    @property
    def unique_score(self):
        AS = self.get_tag('AS', missing=-1)
        XS = self.get_tag('XS', missing=0)
        return AS - XS

    @property
    def edit_distance(self):
        return self.get_tag('NM', missing=9999)

    def __len__(self):
        return self.template_length

    def get_tag(self, *args, missing=None, **kwargs):
        """
        Override `pysam.AlignedSegment.get_tag` method to handle missing tag.
        """
        try:
            return self._pysam_object.get_tag(*args, **kwargs)
        except KeyError:
            return missing

    def get_cigar_operation_length(self, op):
        return dict(self.cigartuples).get('MIDNSHP=XB'.index(op), 0)

    def get_cigar_operation_ratio(self, op):
        cigar_total_length = self.infer_read_length()
        if cigar_total_length is not None:
            length = self.get_cigar_operation_length(op)
            return length / cigar_total_length
        return None

class PysamPileupRead(PysamObjectProxy, PileupReadFilterMixin):
    """
    A proxy class for pysam.PileupRead.
    """
    pass

class PysamPileupColumn(PysamObjectProxy, PileupColumnFilterMixin):
    """
    A proxy class for pysam.PileupColumn.
    """
    pass
