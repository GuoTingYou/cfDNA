from os.path import realpath, abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

import argparse
import gzip

import numpy as np
#import pandas as pd

fractions = np.array([float(f) for f in sys.stdin])
print(fractions.mean(), fractions.std())
print(np.percentile(fractions, [25, 50, 75]))
