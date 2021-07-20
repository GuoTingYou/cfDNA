from ..base import GenomeBase, GenomeMeta, StandardBedSix
from ..base import ReadFileBase, ReadFileMeta

#-------------------------------- Content Class --------------------------------#

class DNaseHotspot(GenomeBase, metaclass = StandardBedSix,
      defaults = (None,)):
    pass

class NarrowPeak(GenomeBase, metaclass = StandardBedSix,
      expand_fields = ('signal_value', 'p_value', 'q_value', 'peak'),
      type_assert = {
          'signal_value': float,
          'p_value': float,
          'q_value': float,
          'peak': int,
      }):
    pass

class BroadPeak(GenomeBase, metaclass = StandardBedSix,
      expand_fields = ('signal_value', 'p_value', 'q_value'),
      type_assert = {
          'signal_value': float,
          'p_value': float,
          'q_value': float,
      }):
    pass

class GappedPeak(GenomeBase, metaclass = StandardBedSix,
      expand_fields = ('thick_start', 'thick_end', 'item_rgb',
                       'block_count', 'block_sizes', 'block_starts',
                       'signal_value', 'p_value', 'q_value'),
      type_assert = {
          'thick_start': int, 'thick_end': int, 'block_count': int,
          'block_sizes':  lambda s: tuple(int(n) for n in s.split(',')),
          'block_starts': lambda s: tuple(int(n) for n in s.split(',')),
          'signal_value': float, 'p_value': float, 'q_value': float,
      }):
    pass

#------------------------------- Read File Class -------------------------------#

class ReadDNaseHotspot(ReadFileBase, metaclass = ReadFileMeta,
                       content_class = DNaseHotspot):
    pass

class ReadNarrowPeak(ReadFileBase, metaclass = ReadFileMeta,
                     content_class = NarrowPeak):
    pass

class ReadBroadPeak(ReadFileBase, metaclass = ReadFileMeta,
                    content_class = BroadPeak):
    pass

class ReadGappedPeak(ReadFileBase, metaclass = ReadFileMeta,
                     content_class = GappedPeak):
    pass
