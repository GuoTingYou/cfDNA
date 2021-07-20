import os

def apply(self, **kwargs):
    """
    This function `apply` is a bounded method of namedtuple derived class,
    so object `self` should have _replace method inherited from namedtuple.

    If `value` is a callable function, then apply this function to the given
    attribute `name` of self. That is, `value(self.name)`.

    If `value` isn't a callable object, then just replace the original value
    of `self.name` with this new given `value`.
    """
    kwargs = {name: value(getattr(self, name)) if callable(value) else value
              for name, value in kwargs.items()}
    return self._replace(**kwargs)

def empty(NT): # NT: namedtuple
    """
    Return an empty instance of given namedtuple derived class `NT`, so
    `NT` should have _make method inherited from namedtuple.
    e.g. `SNP(chrom=None, pos=None, genotype=None)`
    """
    return NT._make([None] * len(NT._fields))

def outfile(inFile, outPath=None, prefix=None, suffix=None):
    """
    Return absolute path of output file. (default: gzip compressed)
    If `outPath` is not given, then use current working directory as `outPath`.
    If `prefix` is not given, then get prefix from the left most "dot" of `inFile`.
    If `suffix` is not given, then use gzip compressed file in default.
    """
    outPath = os.getcwd() if outPath is None else outPath
    prefix = 1 if prefix is None else prefix
    if isinstance(prefix, str):
        prefix = prefix
    elif isinstance(prefix, int):
        prefix = '.'.join(os.path.basename(inFile).split('.')[:prefix])
    else:
        raise ValueError(f'Invalid value "{prefix}" for `prefix`')
    suffix = 'txt' if suffix is None else suffix
    return os.path.join(outPath, f'{prefix}.{suffix}')

def strand_converter(strand):
    forward_sign = ('+', 'forward', 'positive', 'plus')
    reverse_sign = ('-', 'reverse', 'negative', 'minus')

    if not isinstance(strand, (str, type(None))):
        raise TypeError(f'`strand` should be of type `str`, '
                        f'but got {type(strand)!r}.')
    elif strand in (None, 'None', '.'):
        return None
    elif strand.lower() in forward_sign:
        return '+'
    elif strand.lower() in reverse_sign:
        return '-'
    else:
        raise ValueError(f'Valid value for strand should be one of the '
                         f'following: {list(zip(forward_sign, reverse_sign))}')

def pipeline(inputs, *steps):
    index = 0
    task = steps[index]
    chain_map = map(task, inputs)
    for task in steps[1:]:
        chain_map = map(task, chain_map)
    return chain_map

def throw(exception):
    raise exception
