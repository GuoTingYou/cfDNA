#!/home/guotingyou/anaconda3/bin/python
### Filename: baseNS.py ###
# -*- Author: GuoTingYou -*- #
# -*- Coding: utf-8 -*- #

import os
import sys
import time

import re
import gzip
from collections import defaultdict, namedtuple

import numpy as np
from scipy.signal import fftconvolve, find_peaks

### Class and Function ###

class ReadPileup:
    def __init__(self, file):
        self.file = file
        self._open(file)
        self._process()

    def __iter__(self):
        for line in self.content:
            chrom, pos, ref, depth = line.split('\t')[:4]
            yield chrom, int(pos), ref, int(depth)

    def _open(self, file):
        self.content = (line.strip() for line in gzip.open(file, 'rt')
                                     if not line.startswith('#'))

    def _process(self):
        filename = os.path.basename(self.file)
        sample, gene = filename.split('.')[:2]
        self.sample = sample
        self.gene = gene

    def get_depth_array(self, depth_mean=None):
        depth_array = np.array([depth for *_, depth in self])
        depth_mean = depth_mean if depth_mean else depth_array.mean()
        Y = depth_array / depth_mean

        # smoothen
        smooth = 147
        Y = fftconvolve(Y, np.ones(smooth) / smooth, 'valid')
        n = (Y.size - 1) / 2 # TSS flank size
        X = np.linspace(-n, n, Y.size)
        return np.round(X), Y

class ReadNS:
    def __init__(self, file):
        self.file = file
        self._open(file)

    def __iter__(self):
        for line in self.content:
            yield ContentNS(line, self.params)

    def _open(self, file):
        params = [line.strip() for line in gzip.open(file, 'rt')
                               if line.startswith('##')]
        self.params = self._process_params(params)
        self.content = (line.strip() for line in gzip.open(file, 'rt')
                                     if not line.startswith('#'))

    def _process_params(self, params):
        params = [line.strip('## ').split(':') for line in params]
        params = {name: int(value.strip()) for name, value in params}
        return params

class ContentNS:
    def __init__(self, line, params):
        self._process(line, params)

    def _process(self, line, params):
        index, record = line.split('\t')
        ns_array = np.array([ns for ns in record.split('|')], dtype=float)
        self.index = self._process_index(index)
        self.ns_array = self._process_NS_array(ns_array, params)
        self.is_bad_NS = (ns_array.mean() > 4)

    def _process_index(self, index):
        sample_name, chrom, tss, gene_symbol = index.split('::')
        IDX = namedtuple('Index', 'sample, chrom, tss, gene')
        return IDX(sample_name, chrom, int(tss), gene_symbol)

    def _process_NS_array(self, ns_array, params):
        window = params['window']
        smooth = params['smooth']

        # convolution padding border
        border = int(smooth/2) if smooth % 2 == 0 else int((smooth-1)/2)

        # ns x coordinate up- and downstream of TSS,
        # e.g. from -1000 to 1000 if size is 2001.
        coordx = (ns_array.size - 1) / 2

        # align center to 0, which is TSS
        graphx = np.linspace(-coordx, coordx, ns_array.size)
        # x coordinate after convolution (smoothen)
        graphx = graphx[border:-border]

        # smoothen ns by average
        kernel = np.ones(smooth) / smooth
        graphy = fftconvolve(ns_array, kernel, 'valid')

        # setting
        NS = namedtuple('NucleosomeScore', 'X, Y')
        return NS(graphx, graphy)

    def find_nucleosome_center(self):
        X, Y = self.ns_array.X, self.ns_array.Y
        P, _ = find_peaks( Y, distance=147, prominence=0.2, height=1)
        V, _ = find_peaks(-Y, distance=147, prominence=0.2, height=(None, -1))
        self.Peak = X[P]
        self.Valley = X[V]

    def render_score(self, ax, **kwargs):
        X, Y = self.ns_array.X, self.ns_array.Y
        ax.plot(X, Y, **kwargs)

    def render_position(self, ax, y, color, valley=False, **kwargs):
        full = np.full_like
        peakX, peakY = self.Peak, full(self.Peak, y)
        ax.scatter(peakX, peakY, c=color, **kwargs)
        if valley:
            valleyX, valleyY = self.Valley, full(self.Valley, y)
            ax.scatter(valleyX, valleyY, c='', edgecolors=color, **kwargs)
