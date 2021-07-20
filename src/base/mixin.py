import re

from .proxy import PysamAlignedSegment

class ComparableMixin:
    """
    A mixin class to make a genomic locus or region
    comparable between one another.
    """
    def __eq__(self, other):
        if not self.is_empty and not other.is_empty:
            return self._location == other._location
        return False

    def __lt__(self, other):
        if not self.is_empty and not other.is_empty:
            return self._location < other._location
        return False

class AlignedStartMixin(PysamAlignedSegment):
    """
    Find read start site for Nucleosome Score analysis.
    """
    ### Public attributes or methods ###

    def __init__(self, pysam_read):
        super().__init__(pysam_read)
        self.aligned_start = self._get_read_start()

    def relative_start_to(self, reference):
        """
        Adjust read start according to the argument `reference`.
        e.g. A gene's transcription start site.

        The result tell how many base pairs the read start site is
        upstream (negative) or downstream (positive) relative to the `reference`.
        """
        return self.aligned_start - reference

    ### Private attributes or methods ###

    @property
    def _is_cigar_allM(self):
        # consider soft- and hard-clip
        pattern = re.compile(r'^\d+M$')
        if self.cigarstring:
            return True if pattern.match(self.cigarstring) else False
        return False

    def _get_read_start(self):
        """
        Get read start site.
        If all bases of the read is aligned to reference (all Match/Mismatch),
        determined by CIGAR, start is the terminal of the read.
        Otherwise, start is the first aligned base of the read.

        Returns
        -------
        start : [int] read start site
            The position of terminal aligned base.

        Note
        ----
        The coordinate of is 0-based, half-open, end-exclusive.
        To convert the coordinate to that of BAM (1-based, fully-closed),
        start += 1, end = end.
        """
        if self._is_cigar_allM:
            start = self.reference_end if self.is_reverse \
               else self.reference_start + 1
        else:
            align = self.query_alignment_end - 1 if self.is_reverse \
               else self.query_alignment_start # -1 for end point exculsive
            start = self._get_aligned_start(align)
        return start

    def _get_aligned_start(self, align):
        """
        The function is called when any base of the read is not aligned
        to the reference. It return the reference position of the last
        aligned base.

        Parameter
        ---------
        read : [pysam.Alignment object] (required)

        align : [int] first aligned position of the query (required)

        Return
        ------
        rpos : [int, 1-based] reference position of the last aligned base
        """
        for qpos, rpos in self.get_aligned_pairs():
            if rpos and qpos == align:
               return rpos + 1 # to 1-based
        return 0 # something wrong
