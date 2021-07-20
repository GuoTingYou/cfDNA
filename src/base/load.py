#!/home/guotingyou/anaconda3/bin/python
### Filename: read.py ###
# -*- Author: GuoTingYou -*- #
# -*- Coding: utf-8 -*- #

import gzip
import typing

from abc import ABC, ABCMeta
from collections.abc import Iterable
from itertools import takewhile

from .. import util

#---------------------------------- Metaclass ----------------------------------#

class ReadFileMeta(ABCMeta):

    def __new__(cls, clsname, bases, attrs, *, content_class):
        attrs['content_class'] = content_class
        return type(clsname, bases, attrs)

#---------------------------------- Baseclass ----------------------------------#

class ReadFileBase(ABC):

    def __init__(self, file, *args, **kwargs):
        self.file = file
        self.generator = iter(self)

    def __iter__(self):
        if isinstance(self.file, str) or hasattr(self.file, 'read'):
            return util.pipeline(
                self.content,
                lambda args: self.content_class(*args),
                lambda obj: obj.typed)
        elif isinstance(self.file, typing.Generator):
            return self.content
        elif isinstance(self.file, (tuple, list)):
            return iter(self.file)

    def __next__(self):
        try:
            return next(self.generator)
        except StopIteration:
            return util.empty(self.content_class)

    def _open(self):
        if isinstance(self.file, str):
            if self.file.endswith('gz'):
                return gzip.open(self.file, 'rt')
            else:
                return open(self.file, 'r')
        elif hasattr(self.file, 'read'):
            return self.file
        elif isinstance(self.file, typing.Generator):
            return self.file

    def close(self):
        if isinstance(self.generator, typing.Generator):
            self.generator.close()

    @property
    def params(self):
        condition = lambda line: line.startswith('##')
        params = (line for line in takewhile(condition, self._open()))
        return dict(line.strip('#').strip().split(':') for line in params)

    @property
    def fields(self):
        condition = lambda line: line.startswith('#') and not line.startswith('##')
        fields = (line for line in takewhile(condition, self._open()))
        return ' '.join(next(fields).strip('#').strip().split('\t'))

    @property
    def content(self):
        content = self._open()
        if hasattr(content, 'read'):
            return (line.strip().split('\t') for line in content
                                             if not line.startswith('#'))
        else:
            return content
