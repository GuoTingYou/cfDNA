import numpy as np

from ..base import GenomeBase, GenomeMeta
from ..base import ReadFileBase, ReadFileMeta
from ..base import SignalArray

class NucleosomeScoreArray(SignalArray):

    def __new__(cls, nsrecord, **kwargs):
        if isinstance(nsrecord, str):
            ins = super().__new__(cls, np.array(nsrecord.split('|'), dtype=np.float64))
        else:
            ins = super().__new__(cls, nsrecord, **kwargs)
        ins.coverage_too_low = ins.mean() > 4
        return ins

class TssRegion(GenomeBase, metaclass = GenomeMeta,
      expand_fields = ('strand', 'symbol', 'transcription_start_site', 'nucleosome_score'),
      type_assert = {
          'transcription_start_site': int,
          'nucleosome_score': NucleosomeScoreArray,
      }):
    def infer_nucleosome_position(self):
        peaks_up, peaks_down = self.nucleosome_score.find_peak()

        size = self.nucleosome_score.size
        half = int(size / 2)
        # rel for relative
        # ref for reference
        rel_position_x = np.linspace(-half, half, size)
        ref_position_x = np.linspace(self.start, self.end, size)
        # nucleosome position
        rel_nucleosome = rel_position_x[peaks_up]
        ref_nucleosome = ref_position_x[peaks_up]
        # open chromotin region
        rel_ochromotin = rel_position_x[peaks_down]
        ref_ochromotin = ref_position_x[peaks_down]
        # setting
        setattr(self, 'rel_position_x', rel_position_x)
        setattr(self, 'ref_position_x', ref_position_x)
        setattr(self, 'rel_nucleosome', rel_nucleosome)
        setattr(self, 'ref_nucleosome', ref_nucleosome)
        setattr(self, 'rel_ochromotin', rel_ochromotin)
        setattr(self, 'ref_ochromotin', ref_ochromotin)

    def render_nucleosome_score(self, ax, reference=False, **kwargs):
        if reference:
            ax.plot(self.ref_position_x, self.nucleosome_score, **kwargs)
        else:
            ax.plot(self.rel_position_x, self.nucleosome_score, **kwargs)

    def render_nucleosome_position(self, ax, y, color, reference=False, **kwargs):
        kwargs.update(c=color)
        nucleosome_y = np.full_like(self.rel_nucleosome, y)
        if reference:
            ax.scatter(self.ref_nucleosome, nucleosome_y, **kwargs)
        else:
            ax.scatter(self.rel_nucleosome, nucleosome_y, **kwargs)

    def render_ochromotin_position(self, ax, y, color, reference=False, **kwargs):
        kwargs.update(c='', edgecolors=color)
        ochromotin_y = np.full_like(self.rel_ochromotin, y)
        if reference:
            ax.scatter(self.ref_ochromotin, ochromotin_y, **kwargs)
        else:
            ax.scatter(self.rel_ochromotin, ochromotin_y, **kwargs)

class ReadTssRegion(ReadFileBase, metaclass = ReadFileMeta, content_class = TssRegion):
    pass
