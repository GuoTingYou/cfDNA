from ..base import GenomeBase, GenomeMeta
from ..base import ReadFileBase, ReadFileMeta
from ..base import _float

#-------------------------------- Content Class --------------------------------#

class MAF(GenomeBase, metaclass = GenomeMeta,
      expand_fields = ('depth', 'major_allele', 'minor_allele',
                       'fwdDepthA', 'rvsDepthA', 'fwdDepthB', 'rvsDepthB',
                       'major_allele_frequency', 'minor_allele_frequency'),
      type_assert = {
          'depth': int,
          'fwdDepthA': int, 'rvsDepthA': int, # forward and reverse depth of major allele
          'fwdDepthB': int, 'rvsDepthB': int, # forward and reverse depth of minor allele
          'major_allele_frequency': _float,   # allele frequency may be NaN because this
          'minor_allele_frequency': _float,   # locus is indel, homozygous or depth is 0.
      }):
    @property
    def is_indel(self):
        return (self.major_allele in ('+', '-') or
                self.minor_allele in ('+', '-'))

    def get_regionID(self, region_size):
        return int(self.pos / region_size)

#------------------------------- Read File Class -------------------------------#

class ReadMAF(ReadFileBase, metaclass = ReadFileMeta, content_class = MAF):
    pass
