import sys
from collections import namedtuple
from functools import partial
from numbers import Real

from .. import util
from .mixin import ComparableMixin
from .type import Chromosome

#---------------------------------- Baseclass ----------------------------------#

class GenomeBase(ComparableMixin):
    """
    Coordinate: 0-based, half-open, end-exclusive.
    """
    @property
    def size(self):
        return self.end - self.start

    @property
    def is_empty(self):
        return self.chrom is None# or self.chrom.symbol is None
    is_null = is_sentinel = is_empty

    @classmethod
    def sentinel(cls):
        return util.empty(cls)

    def __hash__(self):
        return hash(self._location)

    def __contains__(self, other):
        if self.chrom == other.chrom:
            return self.start <= other.start and other.end <= self.end
        return False

    def signal_from_bigwig(self, bwFile):
        return bwFile.values(self.chrom.symbol, self.start, self.end, numpy=True)

    def overlap(self, other):
        """
        Return:
            result: [float]
                fraction of self region that is overlapped with other.
                If fraction is not `None` and result is not greater than
                `fraction`, return 0.
        """
        if self.chrom != other.chrom:
            return 0
        overlapped = min(self.end, other.end) - max(self.start, other.start)
        overlapped = overlapped / self.size
        return overlapped if overlapped > 0 else 0

    def intersect(self, other, *, fraction=None, reciprocal=False, either=False):
        errmsg_fraction = (
            f"Parameter `fraction` should be a `real number` or (`tuple`, `list`) "
            f"of two real numbers, e.g. (f_self, f_other). Not {type(fraction)!r}."
        )
        if reciprocal and either:
            raise ValueError("Cannot set both `reciprocal` and `either` to True.")

        if fraction is None:
            f_self, f_other = 0.0, 0.0
        elif isinstance(fraction, Real) and reciprocal:
            f_self, f_other = fraction, fraction
        elif isinstance(fraction, Real):
            f_self, f_other = fraction, 0.0
        elif isinstance(fraction, (tuple, list)) and len(fraction) == 2:
            f_self, f_other = fraction
        else:
            raise ValueError(errmsg_fraction)
        return self._intersect(other, (f_self, f_other), either)

    def ison(self, chrom):
        return self.chrom == chrom 

    def fetch_param(self):
        return {
            'contig': self.chrom.symbol,
            'start':  self.start,
            'stop':   self.end,
        }
    def pileup_param(self):
        return {
            'contig':   self.chrom.symbol,
            'start':    self.start,
            'stop':     self.end,
            'truncate': True,
            'stepper':  'all',
        }
    def print(self, *args, fields=None, **kwargs):
        if fields is None:
            args = tuple(getattr(self, attr) for attr in self._fields) + args
        elif isinstance(fields, (tuple, list)):
            args = tuple(getattr(self, attr) for attr in fields) + args
        else:
            raise TypeError('Argument `fields` only accept "tuple" or "list".')
        print(*args, **kwargs)

    @property
    def _location(self):
        return self.chrom._location, self.start, self.end

    def _intersect(self, other, fraction, either):
        f_self, f_other = fraction
        o_self = self.overlap(other)
        o_other = other.overlap(self)
        self_pass = o_self if o_self >= f_self else 0
        other_pass = o_other if o_other >= f_other else 0

        if either and (self_pass or other_pass):
            return o_self, o_other
        return self_pass, other_pass

#---------------------------------- Metaclass ----------------------------------#

class GenomeMeta(type):
    meta_fields = ('chrom', 'start', 'end')
    meta_typed = {'chrom': Chromosome, 'start': int, 'end': int}

    def __new__(cls, clsname, bases, attrs, *,
                defaults=None, expand_fields=None, type_assert=None):
        defaults      = tuple() if defaults is None else defaults
        expand_fields = tuple() if expand_fields is None else expand_fields
        type_assert   = dict()  if type_assert is None else type_assert
        typed = cls.meta_typed.copy(); typed.update(type_assert)
        bases += (namedtuple(clsname, cls.meta_fields + expand_fields,
                             defaults = defaults),)
        attrs['typed'] = property(partial(util.apply, **typed))
        attrs['apply'] = util.apply
        typed.clear()
        return type.__new__(cls, clsname, bases, attrs)
