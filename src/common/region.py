from ..base import GenomeBase
from ..base import GenomeMeta
from ..base import ReadFileBase
from ..base import ReadFileMeta

class RegionalFeature(GenomeBase, metaclass = GenomeMeta,
      expand_fields = ('name', 'depth', 'nloci', 'fraction'),
      type_assert = {
          'depth': int,
          'nloci': int,
          'fraction': float,
      }):
    pass

class RegionalFeatureZ(GenomeBase, metaclass = GenomeMeta,
      expand_fields = ('name', 'depth', 'nloci', 'fraction', 'zscore'),
      type_assert = {'name': int, 'depth': int, 'nloci': int,
                     'fraction': float, 'zscore': float}):
    pass

class ReadRegionalFeature(ReadFileBase, metaclass = ReadFileMeta,
                          content_class = RegionalFeature):
    pass

class ReadRegionalFeatureZ(ReadFileBase, metaclass = ReadFileMeta,
                           content_class = RegionalFeatureZ):
    pass
