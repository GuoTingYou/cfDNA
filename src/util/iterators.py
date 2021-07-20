from abc import abstractmethod
from collections.abc import Iterable, Iterator

class DoubleIterator(Iterable):
    """
    Descripsion
    -----------
    The class is an abstract class for iterating two sorted file simultaneously
    and yielding item in sorted order. The iterating blueprint is given in its
    `__iter__` method, and detailed yielding rule should be defined in its
    derived class.

    Parameters
    ----------
    `iter1` and `iter2` should be instance of iterator which contains both
    `__iter__` and `__next__` method. Additionally, both `item1` and `item2`
    from `iter1` and `iter2` respectively should be comparable. That is,
    metaclass of themselves or their baseclass should be `GenomeComparableMeta`
    in ../base/genome.py module.

    Examples
    --------
    File1: Locus(chrom="chr1", pos="100")
           Locus(chrom="chr2", pos="100")
           Locus(chrom="chrX", pos="500")

    File2: Locus(chrom="chrM", pos="200")
           Locus(chrom="chr1", pos="200")
           Locus(chrom="chr1", pos="300")

    Yield: Locus(chrom="chrM", pos="200")
        -> Locus(chrom="chr1", pos="100")
        -> Locus(chrom="chr1", pos="200")
        -> Locus(chrom="chr1", pos="300")
        -> Locus(chrom="chr2", pos="100")
        -> Locus(chrom="chrX", pos="500")
    """
    def __init__(self, iter1, iter2):
        self.iter1 = self._check_iterable(iter1)
        self.iter2 = self._check_iterable(iter2)
        self.item1 = next(self.iter1)
        self.item2 = next(self.iter2)

    def __iter__(self):
        while not self.item1.is_empty \
          and not self.item2.is_empty:
            yield self.compare()

        while not self.item1.is_empty:
            yield self.case_iter1_unfinished(self.item1)
            self.item1 = next(self.iter1)

        while not self.item2.is_empty:
            yield self.case_iter2_unfinished(self.item2)
            self.item2 = next(self.iter2)

    def compare(self):
        if self.item1 == self.item2:
            result = self.case_equal(self.item1, self.item2)
            self.item1 = next(self.iter1)
            self.item2 = next(self.iter2)
            return result

        elif self.item1 < self.item2:
            result = self.case_smaller(self.item1, self.item2)
            self.item1 = next(self.iter1)
            return result

        elif self.item1 > self.item2:
            result = self.case_larger(self.item1, self.item2)
            self.item2 = next(self.iter2)
            return result

        else:
            raise TypeError(
                f'{self.item1} of type {type(self.item1)} or '
                f'{self.item2} of type {type(self.item2)} is '
                 'uncomparable.')

    @abstractmethod
    def case_equal(self, item1, item2):
        pass
    @abstractmethod
    def case_smaller(self, item1, item2):
        pass
    @abstractmethod
    def case_larger(self, item1, item2):
        pass
    @abstractmethod
    def case_iter1_unfinished(self):
        pass
    @abstractmethod
    def case_iter2_unfinished(self):
        pass

    def _check_iterable(self, iter_obj):
        if isinstance(iter_obj, Iterator):
            return iter_obj
        raise TypeError(
            f'Input `{iter_obj}` should be instance of iterator '
             'which contains both `__iter__` and `__next__` method.')
