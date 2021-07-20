from os.path import dirname, join
from numpy import nan

from .mixin import ComparableMixin


global GENOME_SIZE, CANONICAL_CHROMOSOMES

def _load_genome_size(filename):
    with open(filename) as f:
        f = (line.strip().split('\t') for line in f)
        return {chrom: int(size) for chrom, size in f}

GENOME_SIZE = _load_genome_size(join(dirname(__file__), '../../cfg/genome.hg19'))
CANONICAL_CHROMOSOMES = tuple(f'chr{N}' for N in ['M', *range(1, 23), 'X', 'Y'])

class Chromosome(ComparableMixin):

    def __init__(self, symbol):
        self.symbol = symbol
        self.min_pos = 0
        self.max_pos = GENOME_SIZE.get(self.symbol, -1)
        self.is_empty = symbol is None

    def __getattr__(self, attr):
        return getattr(self.symbol, attr)

    def __len__(self):
        return GENOME_SIZE.get(self.symbol, -1)

    def __repr__(self):
        return f'{type(self).__name__}({self.symbol}:{self.min_pos}-{self.max_pos})'

    def __str__(self):
        return self.symbol

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.symbol == other.symbol
        return False

    def __le__(self, other):
        return self == other or self < other

    def __hash__(self):
        return hash(self.symbol)

    @property
    def is_canonical(self):
        return self.symbol in CANONICAL_CHROMOSOMES

    @property
    def _location(self):
        try:
            return (tuple(GENOME_SIZE.keys()).index(self.symbol),)
        except ValueError:
            return self._missing()

    def _missing(self):
        return (len(GENOME_SIZE),) + tuple(ord(char) for char in self.symbol)

def _float(data):
    """
    Allele frequency may be NaN because this
    locus is indel, homozygous or depth is 0.
    Corresponding to `cal_allele_frequencies`
    method of class `Pileup` in
    `../pileup/pileup.py module`.
    Thus, convert it to numpy.nan for
    simplifying computation
    """
    try:
        return float(data)
    except ValueError as err:
        if data in ('None', 'NA', 'nan'):
            return nan
        else:
            raise ValueError(err)
