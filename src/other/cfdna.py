#!/home/guotingyou/anaconda3/bin/python
### Filename: cfdna.py ###
# -*- Author: GuoTingYou -*- #
# -*- Coding: utf-8 -*- #

import pysam

class ComparableMixin:
    """
    A mixin class for PysamAlignedSegmentProxy to make pysam.AlignedSegment
    objects comparable between one another.
    """
    contig = dict((chrom, n) for n, chrom in enumerate(
                   f'chr{N}' for N in ['M', *range(1, 23), 'X', 'Y']))
    @property
    def _chrom_start_end(self):
        return self.contig.get(self.chrom, 25), self.start, self.end

    def __eq__(self, other):
        return self._chrom_start_end == other._chrom_start_end

    def __lt__(self, other):
        return self._chrom_start_end < other._chrom_start_end

class FilterMixin:
    """
    A mixin class for PysamAlignedSegmentProxy to apply some basic filterations
    on pysam.AlignedSegment objects comparable between one another.
    """
    def apply_filters(self, *, MQ = 30,
                      exclude_supplementary = True,
                      exclude_duplicate = True,
                      exclude_secondary = True,
                      exclude_qcfail = True):
        return (not self.is_sup(exclude_supplementary) and
                not self.is_dup(exclude_duplicate) and
                not self.is_sec(exclude_secondary) and
                not self.is_qcf(exclude_qcfail) and
                self.mapping_quality >= MQ)

    def is_sup(self, exclude_supplementary):
        return self.is_supplementary if exclude_supplementary else False

    def is_dup(self, exclude_duplicate):
        return self.is_duplicate if exclude_duplicate else False

    def is_sec(self, exclude_secondary):
        return self.is_secondary if exclude_secondary else False

    def is_qcf(self, exclude_qcfail):
        return self.is_qcfail if exclude_qcfail else False

class PysamAlignedSegmentProxy(ComparableMixin, FilterMixin):
    """
    A proxy class of pysam.AlignedSegment.
    """
    def __init__(self, pysam_read):
        self.pread = pysam_read
        self.chrom = self.reference_name
        self.start = self.reference_start
        self.end   = self.reference_end
        self.length = self.query_alignment_length
        self.strand = '-' if self.is_reverse else '+'

    def __str__(self):
        return self.pread.__str__()
    __repr__ = __str__

    def __len__(self):
        return self.template_length

    def __getattr__(self, attr):
        try:
            return getattr(self.pread, attr)
        except AttributeError:
            raise AttributeError(f'{self.__class__.__name__} or pysam.AlignedSegment '
                                 f'object has no attribute {attr}.')

class ReadStartSite(PysamAlignedSegmentProxy):
    """
    Find read start site for Nucleosome Score analysis.
    """
    def __init__(self, pysam_read):
        super().__init__(pysam_read)
        self.aligned_start = self._get_read_start()

    @property
    def is_cigar_allM(self):
        # consider soft- and hard-clip
        pattern = re.compile(r'^\d+M$')
        cigar = self.cigarstring
        if cigar is None: return False
        return True if pattern.match(cigar) else False

    def relative_start_to(self, pos):
        """
        Adjust read start according to the argument `pos`.
        e.g. A gene's transcription start site.

        The result tell how many base pairs the read start site is
        upstream (negative) or downstream (positive) relative to the `pos`.
        """
        return self.start - pos

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
        if self.is_cigar_allM:
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

#    bamFile = pysam.AlignmentFile(sys.argv[1])
#    for read in map(CellFreeDNA, bamFile):
#        print(read.chrom, read.length, read.strand)
#        print(read)
#        break

# reads = (read for read in map(CellFreeDNA, bamFile) if read.apply_filters(MQ=30))
