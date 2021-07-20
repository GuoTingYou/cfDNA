from ..base import GenomeBase
from ..base import GenomeMeta
from ..base import ReadFileBase
from ..base import ReadFileMeta
from ..base import StandardBedSix
from ..base import StandardBedTwelve
from ..util import pipeline, throw

class GeneBase(GenomeBase):
    @property
    def transcription_start_site(self):
        return self.start if self.strand == '+' \
          else self.end if self.strand == '-' \
          else throw(ValueError(f'Cannot find transcription start site for '
                                f'gene: "{self}".'))
    @property
    def transcription_end_site(self):
        return self.end if self.strand == '+' \
          else self.start if self.strand == '-' \
          else throw(ValueError(f'Cannot find transcription end site for '
                                f'gene: "{self}".'))

    def tss_region(self, extend=1000):
        tss = self.transcription_start_site
        return self.chrom, tss - extend, tss + extend

class GeneText(GeneBase, metaclass = GenomeMeta,
      expand_fields = ('strand', 'symbol'),
      type_assert = {
          'strand': (lambda s: '+' if s.lower() in ('+', 'forward', 'positive')
                          else '-' if s.lower() in ('-', 'reverse', 'negative')
                          else throw(ValueError(f'Invalid strand "{s}"')))
      }):
    pass

class GeneBed6(GeneBase, metaclass = StandardBedSix):
    pass

class GeneBed12(GeneBase, metaclass = StandardBedTwelve):
    pass

class GeneTable(ReadFileBase, metaclass = ReadFileMeta, content_class = GeneText):
    def __init__(self, *args, fmt='txt', **kwargs):
        if fmt == 'txt':
            pass
        elif fmt == 'bed6':
            self.content_class = GeneBed6
        elif fmt == 'bed12':
            self.content_class = GeneBed12
        else:
            raise ValueError(f'Invalid format "{fmt}" for {type(self).__name__}')
        super().__init__(*args, **kwargs)

    def __iter__(self):
        step1 = lambda ln: ln[:len(self.content_class._fields)]
        step2 = lambda ln: [ln[0].split('_')[0], *ln[1:]] # split chrN_globxxx
        step3 = self.content_class._make
        step4 = lambda gene: gene.typed
        return pipeline(self.content, step1, step2, step3, step4)
