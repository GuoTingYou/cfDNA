from ..util import strand_converter

from .genome import GenomeBase
from .genome import GenomeMeta
from .load import ReadFileBase
from .load import ReadFileMeta

BED_FIELD_NAMES = (
    'chrom',
    'start',
    'end',
    'name',
    'score',
    'strand',
    'thick_start',
    'thick_end',
    'item_rgb',
    'block_count',
    'block_sizes',
    'block_starts',
)

class StandardBedSix(GenomeMeta):
    def __new__(cls, *args, **kwargs):
        cls.meta_fields = tuple(super().meta_fields) + BED_FIELD_NAMES[3:6]
        cls.meta_typed = dict(super().meta_typed.items())
        cls.meta_typed.update({
            BED_FIELD_NAMES[4]: float,
            BED_FIELD_NAMES[5]: strand_converter,
        })
        return super().__new__(cls, *args, **kwargs)

class StandardBedTwelve(GenomeMeta):
    def __new__(cls, *args, **kwargs):
        cls.meta_fields = tuple(super().meta_fields) + BED_FIELD_NAMES[3:12]
        cls.meta_typed = dict(super().meta_typed.items())
        cls.meta_typed.update({
            BED_FIELD_NAMES[4]: float,
            BED_FIELD_NAMES[5]: strand_converter,
            BED_FIELD_NAMES[6]: int,
            BED_FIELD_NAMES[7]: int,
            BED_FIELD_NAMES[9]: int,
            BED_FIELD_NAMES[10]: lambda x: [int(y) for y in x.strip(',').split(',')],
            BED_FIELD_NAMES[11]: lambda x: [int(y) for y in x.strip(',').split(',')],
        })
        return super().__new__(cls, *args, **kwargs)

class BedThreeRegion(GenomeBase, metaclass=GenomeMeta):
    pass

class ReadBed(ReadFileBase, metaclass=ReadFileMeta, content_class=BedThreeRegion):
    pass
