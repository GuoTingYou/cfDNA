from ..base import GenomeBase, GenomeMeta
from ..base import ReadFileBase, ReadFileMeta
from ..base import _float
from ..util import empty

#-------------------------------- Content Class --------------------------------#

class SNP(GenomeBase, metaclass = GenomeMeta,
      expand_fields = ('genotype', 'allele_A', 'allele_B',
                       'fwdDepthA', 'rvsDepthA', 'fwdDepthB', 'rvsDepthB',
                       'fraction'),
      type_assert = {
          'fwdDepthA': int, 'rvsDepthA': int, # forward and reverse depth of allele A
          'fwdDepthB': int, 'rvsDepthB': int, # forward and reverse depth of allele B
          'fraction': _float,                 # fraction of a single locus
      }):
    def estimate_fraction(self, DepthA, DepthB):
        try:
            return {
                'AA-ab': lambda a, b: b * 2 / (a + b),
                'AB-aa': lambda a, b: (a-b) / (a + b),
                'AA-bb': lambda a, b: b     / (a + b),
            }[self.genotype.split(':')[0]](DepthA, DepthB)
        except ZeroDivisionError:
            return None

    @property
    def is_indel(self):
        return False

    @property
    def background_allele(self):
        return {
            'AA-ab': 'a',
            'AB-aa': '1b',
            'AA-bb': '2a',
        }[self.genotype.split(':')[0]]

    @property
    def interested_allele(self):
        return {
            'AA-ab': '1b',
            'AB-aa': 'a',
            'AA-bb': '2b',
        }[self.genotype.split(':')[0]]

    @property
    def depth(self):
        return self.depth_sum()

    def depth_of(self, allele):
        if allele == 'a':
            return self.depth_sum(which='deptha') - self.depth_sum(which='depthb')
        elif allele == '1b':
            return self.depth_sum(which='depthb') * 2
        elif allele == '2b':
            return self.depth_sum(which='depthb')
        elif allele == '2a':
            return self.depth_sum(which='deptha')
        else:
            raise ValueError(f'Invalid argument {allele} for allele. Usage: '
                             'self.depth_of(allele=self.interested_allele) or '
                             'self.depth_of(allele=self.background_allele)')

    def check_fraction(self, *, fmin=0, fmax=1):
        return fmin < self.fraction < fmax

    def check_depth(self, *, at_least):
        depthA = self.depth_sum('A')
        depthB = self.depth_sum('B')
        less = min(depthA, depthB)
        more = max(depthA, depthB)
        return (less > more / 50) and (less >= at_least)

    def get_regionID(self, region_size):
        return int(self.pos / region_size)

    def depth_sum(self, which='depth'):
        which = which if which in ('fwd', 'rvs', 'deptha', 'depthb', 'depth') else 'depth'
        return sum([getattr(self, a) for a in self._fields if which in a.lower()])

#------------------------------- Read File Class -------------------------------#

class ReadSNP(ReadFileBase, metaclass = ReadFileMeta, content_class = SNP):
    pass
