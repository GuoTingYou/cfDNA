from ..util import DoubleIterator

class IntersectFeatureSelector(DoubleIterator):
    def __init__(self, iter1, iter2, *,
                 fraction=None, reciprocal=False, either=False):
        super().__init__(iter1, iter2)
        self.intersect_kw = {
            'fraction': fraction,
            'reciprocal': reciprocal,
            'either': either,
        }

    def __iter__(self):
        for features in super().__iter__():
            if features is None:
                continue
            elif isinstance(features, tuple):
                yield features
            else:
                yield from features

    def case_equal(self, feature1, feature2):
        return feature1, feature2, 1, 1

    def case_smaller(self, feature1, feature2):
        o1, o2 = feature1.intersect(feature2, **self.intersect_kw)
        while all((o1, o2)): # no 0 in (o1, o2)
            yield feature1, feature2, o1, o2
            feature2 = self.item2 = next(self.iter2)
            o1, o2 = feature1.intersect(feature2, **self.intersect_kw)

    def case_larger(self, feature1, feature2):
        o1, o2 = feature1.intersect(feature2, **self.intersect_kw)
        if all((o1, o2)): # no 0 in (o1, o2)
            return feature1, feature2, o1, o2

    def case_iter1_unfinished(self, feature1):
        return feature1, self.item2, 0, 0

    def case_iter2_unfinished(self, feature2):
        self.iter2.close()
        return
