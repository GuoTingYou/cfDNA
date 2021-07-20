#!/home/guotingyou/anaconda3/bin/python
### Filename: graph_NS.py ###
# -*- Author: GuoTingYou -*- #
# -*- Coding: utf-8 -*- #

import os
import sys
import time

import re
import gzip
import argparse
from functools import wraps
from itertools import takewhile
from collections import defaultdict, namedtuple

import numpy as np
from scipy.signal import fftconvolve, find_peaks

import matplotlib as mpl; mpl.use("PDF")
import matplotlib.pyplot as plt
import seaborn as sns

youyou = "/zfssz2/ST_MCHRI/BIGDATA/USER/guotingyou/youyou/"
tortoise = youyou + "fetal_fraction/"
sys.path.append(tortoise)
import tortoise as tt

### I/O Interface and Usage ###

def argParse(N):
    script = sys.argv[0]
    usage = "Not enough arguments were given. Expected {1}\n{0}".format(
            " ".join(["Usage: python", script]), N)
    description = "{0:^81s}\n{1:^81s}\n{2:^81s}\n{3:^81s}".format(
                  "Filename: {}".format(script),
                  "-*- Author: GuoTingYou -*-",
                  "-*- Version: youyou.tortoise 0.0.0 -*-",
                  "-*- Coding: utf-8 -*-")
    parser = argparse.ArgumentParser(
             prog = script, usage = usage, description = description,
             formatter_class = argparse.RawDescriptionHelpFormatter)
    help_info = {
        "i": "Input: Nucleosome scores <*.ns.gz> file.",
        "O": "Output: Absolute output path.",
        "o": "Output: Prefix of output name and in TITLE.",
        "t": "Title for figure",
    }
    arguments = {
        "v": dict(action = "version", version = "youyou.tortoise 0.0.0"),
        "i": dict(type=str, help=help_info["i"], nargs="?", required=True),
        "O": dict(type=str, help=help_info["O"], nargs="?", default=os.getcwd()),
        "o": dict(type=str, help=help_info["o"], nargs="?", default=None),
        "t": dict(type=str, help=help_info["t"], nargs="?", default="Nucleosome Score of {}"),
    }
    parser.add_argument("-v", "--version",  **arguments["v"])
    parser.add_argument("-i", "--input",    **arguments["i"])
    parser.add_argument("-O", "--out_path", **arguments["O"])
    parser.add_argument("-o", "--out_name", **arguments["o"])
    parser.add_argument("-t", "--title",    **arguments["t"])
    if len(sys.argv) <= N: parser.print_help(); raise IndexError(usage)
    return parser.parse_args()
try:
    ioargs = argParse(1)
    i_file = ioargs.input # required
    o_path = ioargs.out_path
    o_name = ioargs.out_name
    o_fig  = 'Nucleosome_Score.{}.pdf' if o_name is None \
        else '{}.Nucleosome_Score.{}.pdf'.format(o_name, '{}')
except Exception as err:
    raise Exception(err)
else:
    o_fig = os.path.join(o_path.rstrip('/'), o_fig)
    TITLE = ioargs.title
    WINDOW = 'to_be_update'
    SMOOTH = 'to_be_update'
    MERGE  = 'to_be_update'

### Class and Function ###

def timetask(func):
    @wraps(func)
    def decorator(*args, **kwargs):
        t0 = time.time()
        print("{0} [{1}]".format(time.strftime("%X", time.localtime()), func.__name__))
        result = func(*args, **kwargs)
        print("{0} [{1}]".format(time.strftime("%X", time.localtime()), func.__name__))
        t1 = time.time() - t0
        print("Finished [{1}] consumed: {0} Secs.".format(round(t1, 2), func.__name__))
        return result
    return decorator

class ReadNS:
    def __init__(self, file):
        self.file = file
        self._open(file)

    def __iter__(self):
        for line in self.content:
            yield ContentNS(line, self.params)

    def _open(self, file):
        params = [line.strip() for line in gzip.open(file, 'rt') if line.startswith('##')]
        self.params = self._process_params(params)
        self.content = (line.strip() for line in gzip.open(file, 'rt') if not line.startswith('#'))

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
        # normalize y from 0 to 1
        #graphy = (graphy - min(graphy)) / (max(graphy) - min(graphy))

        # find peaks and valley
        P, _ = find_peaks( graphy, distance=147, prominence=0.2)
        V, _ = find_peaks(-graphy, distance=147, prominence=0.2)
        peak = graphx[P]; valley = graphx[V]

        # setting
        NS = namedtuple('NucleosomeScore', 'X, Y, Peak, Valley')
        return NS(graphx, graphy, peak, valley)

    def render_score(self, ax, **kwargs):
        X, Y = self.ns_array.X, self.ns_array.Y
        ax.plot(X, Y, **kwargs)

    def render_position(self, ax, y, color, **kwargs):
        peakX, valleyX = self.ns_array.Peak, self.ns_array.Valley
        peakY, valleyY = np.full_like(peakX, y), np.full_like(valleyX, y)
        ax.scatter(peakX, peakY, c=color, **kwargs)
        ax.scatter(valleyX, valleyY, c='', edgecolors=color, **kwargs)

def figure_setting(ax1, ax2, N, *, title):
    ax1.spines['left'].set_position(('data', 0))
    ax1.set_xlim(-WINDOW, WINDOW)
    anchor = -.125 if N <= 8 else -.175
    ax1.legend(
        loc = 'lower center',
        ncol = 8,
        bbox_to_anchor = (0, anchor, 1., .102),
        markerscale = 1.2,
        frameon = False,
        framealpha = .5,
        fontsize = 16,
    )
    ax1 = tt.setTitle(ax1, t = title)
    ax1 = tt.setSpine(ax1)
    ax1 = tt.setTicks(ax1, gridv = (0, 1))

    ax2.set(
        xlim = (-WINDOW, WINDOW),
        ylim = (.5, N+.5),
        yticks = [],
    )
    ax2 = tt.setTitle(ax2, x = 'Nucleosome Position (Transcription Start Site: 0)')
    ax2 = tt.setSpine(ax2, visible = (0, 0, 0, 1))
    ax2 = tt.setTicks(ax2, gridv = (1, 0), tickv = (0, 0, 0, 1))
    return ax1, ax2

def main_individual():
    global o_name, WINDOW, SMOOTH
    NSFile = ReadNS(i_file)
    WINDOW = NSFile.params['window']
    SMOOTH = NSFile.params['smooth']
    filtered_file = (ns_record for ns_record in NSFile if not ns_record.is_bad_NS)

    gs_prop = {
        'hspace': 0.25,
        'width_ratios':  [1],
        'height_ratios': [3, 1],
    }
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(25, 15), gridspec_kw=gs_prop)

    data = []
    colors = sns.color_palette('Set1')
    for y, ns_record in enumerate(filtered_file, 1):
        ns_record.render_score(ax1, color=colors[y-1], linewidth=2.5,
                               label=ns_record.index.gene)
        ns_record.render_position(ax2, y, colors[y-1], s=100)
        data.append(ns_record)

    outfig = o_fig.format(ns_record.index.sample)
    title = f'Nucleosome Score of {ns_record.index.sample}'
    ax1, ax2 = figure_setting(ax1, ax2, len(data), title=title)
    fig.savefig(outfig, dpi = 800, transparent = True, bbox_inches = 'tight')

def main_list():
    global o_name, WINDOW, SMOOTH

    _N = 0
    data = defaultdict(list)
    for inFile in open(i_file, 'r'):
        _N += 1
        NSFile = ReadNS(inFile.strip())
        for ns_record in NSFile:
            if not ns_record.is_bad_NS:
                data[ns_record.index.gene].append(ns_record)
    WINDOW = NSFile.params['window']
    SMOOTH = NSFile.params['smooth']

    gs_prop = {
        'hspace': 0.25,
        'width_ratios':  [1],
        'height_ratios': [3, 1],
    }
    #colors = sns.color_palette('Set1')
    RED, YELLOW, GREEN, BLUE, PURPLE = '#d3121b', '#fdb813', '#00a86b', '#1867bb', '#8f00ff'
    colors = [RED, YELLOW, GREEN, BLUE, PURPLE]
    #colors = sns.light_palette(BLUE, _N+1)
    for geneSymbol in data.keys():
        data_individual = data[geneSymbol]
        if len(data_individual) == _N:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(25, 15), gridspec_kw=gs_prop)
            for y, ns_record in enumerate(data_individual, 1):
                c = colors[y-1]
                #get_gene(ns_record)
                ns_record.render_score(ax1, color=c, linewidth=2.5,
                                       label=ns_record.index.sample)
                ns_record.render_position(ax2, y, c, s=100)
            outfig = o_fig.format(ns_record.index.gene)
            title = f'Nucleosome Score of {ns_record.index.gene}'
            ax1, ax2 = figure_setting(ax1, ax2, _N, title=title)
            fig.savefig(outfig, dpi = 800, transparent = True, bbox_inches = 'tight')

def get_gene(ns_record):
    if ns_record.index.gene == 'RHOU':
        day = ns_record.index.sample.split('_')[1]
        X = ns_record.ns_array.X
        Y = ns_record.ns_array.Y
        of = open(f'array.ns.LXQ.{day}.txt', 'w')
        for x, y in zip(X, Y):
            print(x, y, sep='\t', file=of)

### Main ###

if __name__ == "__main__":

    print(os.path.realpath(__file__))
    print("\n{:-^81s}".format(
          " Task Start At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))

    #main_individual()
    main_list()

    print("{:-^81s}\n".format(
          " Task End At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))
    print("Practice_makes_perfect")
