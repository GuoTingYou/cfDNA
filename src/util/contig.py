class Contig:
    # {'chrM': 0, 'chr1': 1, ...}
    CHR_to_NUM = dict((chrom, n) for n, chrom in enumerate(
                       f'chr{N}' for N in ['M', *range(1, 23), 'X', 'Y']))
    # {0: 'chrM', 1: 'chr1', ...}
    NUM_to_CHR = dict((n, chrom) for chrom, n in CHR_to_NUM.items())

    def __init__(self, chrom):
        self._chrom = chrom

    def chrom(self):
        if self._chrom.startswith('chr'):
            return self._chrom
        else:
            try:
                return self.NUM_to_CHR[int(self._chrom)]
            except KeyError:
                raise KeyError(f'Unrecognized chromosome number "{self._chrom}"'
                               f'for `NUM_to_CHR`. Valid number to chromosome'
                               f'mappings are: {NUM_to_CHR}.')

    def number(self):
        if self._chrom.startswith('chr'):
            try:
                return self.CHR_to_NUM[self._chrom]
            except KeyError:
                return self._chrom
