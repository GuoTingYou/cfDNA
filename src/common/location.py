from abc import abstractmethod

from ..base import ReadFileBase

class Location:
    def __init__(self, locFile):
        self.locFile = locFile

    def __iter__(self):
        if isinstance(self.locFile, ReadFileBase):
            return map(self.get_param(), self.locFile)
        else:
            return map(self.get_param_from_bed(), ReadBED(self.locFile))

    @abstractmethod
    def get_param(self):
        pass
    @abstractmethod
    def get_param_from_bed(self):
        pass

class FetchLocation(Location):
    def get_param(self):
        return lambda loc: (loc.fetch_param(), loc)

    def get_param_from_bed(self):
        def _fetch_param(self, location):
            chrom, start, end = location
            param = {
                'contig':   chrom,
                'start':    start-1,
                'stop':     end,
            }
            return param, location
        return _fetch_param

class PileupLocation(Location):
    def get_param(self):
        return lambda loc: (loc.pileup_param(), loc)

    def get_param_from_bed(self):
        def _pileup_param(self, location):
            chrom, start, end = location
            param = {
                'contig':   chrom,
                'start':    start-1,
                'stop':     end,
                'truncate': True,
                'stepper':  'all',
            }
            return param, location
        return _pileup_param

class ReadBED(ReadFileBase):

    def __iter__(self):
        try:
            for line in self.content:
                chrom, start, end, *_ = line.strip().split('\t')
                yield chrom, int(start), int(end)
        except ValueError:
            for line in self.content:
                chrom, pos, *_ = line.strip().split('\t')
                yield chrom, int(pos), int(pos)
