from collections import Counter

import numpy as np
from pandas import Series

from ..base import PysamPileupColumn, GenomeBase

class Pileup(PysamPileupColumn, GenomeBase):

    def __init__(self, pileupcolumn, inclusive):
        super().__init__(pileupcolumn)
        self.inclusive = inclusive
        self.depth = inclusive.count(True)
        self.chrom = pileupcolumn.reference_name
        self.start = pileupcolumn.reference_pos
        self.end = pileupcolumn.reference_pos + 1

    @property
    def pass_reads(self):
        return Series(self.pileups)[self.inclusive]

    @property
    def sequences(self):
        # To reproduce samtools mpileup format, set all of
        # mark_matches, mark_ends and add_indels to True of
        # method `get_query_sequences`.
        return Series(self.get_query_sequences())[self.inclusive]

    @property
    def base_qualities(self):
        return Series(self.get_query_qualities())[self.inclusive]

    def get_share_special_alignments(self, alleles: tuple):
        share, special = alleles
        share_aligns, special_aligns = [], []
        for base, pileupread in zip(self.sequences, self.pass_reads):
            if base == share:
                share_aligns.append(pileupread.alignment)
            elif base == special:
                special_aligns.append(pileupread.alignment)
        return share_aligns, special_aligns

    def depth_of(self, allele, strand=None):
        sequences = self.sequences.to_list()
        forward = sequences.count(allele.upper())
        reverse = sequences.count(allele.lower())
        return (forward if strand in ('+', 'forward')
          else  reverse if strand in ('-', 'reverse')
          else  forward, reverse)

    def insertsizes(self, fmt=str):
        get_insertsize = lambda read: fmt(abs(read.alignment.template_length))
        return self.pass_reads.apply(get_insertsize)

    def insertsizes_of(self, allele, strand=None):
        if strand in ('+', 'forward'):
            index = (self.sequences == allele.upper())
        elif strand in ('-', 'reverse'):
            index = (self.sequences == allele.lower())
        else:
            index = (self.sequences.apply(str.upper) == allele.upper())
        return self.insertsizes()[index].values

    def cal_allele_depths(self, bq_weighted=False):
        counts = Counter(self.sequences.apply(str.upper)).most_common(2)
        counts.extend([('None', 0), ('None', 0)])
        (major_allele, major_depth), (minor_allele, minor_depth), *_ = counts

        if bq_weighted:
            bcar = self._base_call_accuracy_rate()
            major_depth = bcar[self.sequences == major_allele].sum()
            minor_depth = bcar[self.sequences == minor_allele].sum()

        setattr(self, 'major_allele', major_allele)
        setattr(self, 'major_depth',  major_depth)
        setattr(self, 'minor_allele', minor_allele)
        setattr(self, 'minor_depth',  minor_depth)

    def cal_allele_frequencies(self, bq_weighted=False):
        self.cal_allele_depths(bq_weighted=bq_weighted)
        if self.depth == 0:
            MAF, maf = 'NA', 'NA'
        else:
            MAF = self.major_depth / self.depth
            maf = self.minor_depth / self.depth
        setattr(self, 'major_allele_frequency', MAF)
        setattr(self, 'minor_allele_frequency', maf)

    def _base_call_accuracy_rate(self):
        baqs = np.array(self.base_qualities)
        return 1 - np.float_power(10, -baqs / 10)
