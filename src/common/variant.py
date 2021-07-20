import vcf

from ..base import GenomeBase, GenomeMeta, ReadFileBase
from ..util import empty

class Variant(GenomeBase, metaclass = GenomeMeta,
      expand_fields = ('genotype', 'is_pass', 'is_snp', 'is_hetero', 'gt_bases')):

    def is_biallelic(self, other):
        return len(set(self.gt_bases.split('/')) | set(other.gt_bases.split('/'))) == 2

    def is_AAab(self, other):
        return self.is_biallelic(other) and (not self.is_hetero and other.is_hetero)

    def is_ABaa(self, other):
        return self.is_biallelic(other) and (self.is_hetero and not other.is_hetero)

    def is_AAbb(self, other):
        return self.is_biallelic(other) and (not self.is_hetero and not other.is_hetero)

class VCF(vcf.Reader):

    def __init__(self, filename):
        super().__init__(filename=filename)

    def _fetch(self, chrom, start, end):
        for record in super().fetch(chrom, start, end):
            call = record.samples[0]
            flag = False if record.FILTER is None else len(record.FILTER) == 0
            if flag and call.is_het:
                flag = all(depth >= 10 for depth in call.data.AD)
            yield Variant(record.CHROM, record.POS-1, record.end, call["GT"], flag,
                          record.is_snp, call.is_het, call.gt_bases).typed

    def fetch(self, chrom, start=None, end=None):
        records = self._fetch(chrom, start, end)
        while True:
            try:
                record = next(records)
                yield record
            except StopIteration:
                yield empty(Variant)
