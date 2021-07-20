from os.path import realpath, abspath, dirname, join

import os
import sys
sys.path.insert(0, abspath(join(dirname(realpath(__file__)), '..')))

import argparse
import gzip

from src import (
    # class
    ReadRegion,
    ReadNarrowPeak,
    PassSnpSelector,
    # function
    argParse_template,
    timemain,
    outfile,
)
