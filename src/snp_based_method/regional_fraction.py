import os
import sys
import time

from base import (
    ReadFileMeta, LocusMeta, RegionMeta,
    ReadFileBase, GenomeComparableBase
)
from util import empty

#-------------------------------- Content Class --------------------------------#

class RegionalFraction(metaclass=RegionMeta, expand_fields=['nsnps', 'fraction'],
    type_assert = {'nsnps': int, 'fraction': float}):
    pass

#------------------------------- Read File Class -------------------------------#

class ReadRegionalFraction(ReadFileBase,
                           metaclass = ReadFileMeta,
                           content_class = RegionalFraction):
    pass

### Main ###

if __name__ == "__main__":

    print(os.path.realpath(__file__))
    print("\n{:-^81s}".format(
          " Task Start At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))

    print("\n{:^81s}\n".format("~*~ TORTOISE ~*~"))

    print("{:-^81s}\n".format(
          " Task End At: " + time.strftime("%Y-%m-%d %H:%M:%S ", time.localtime())
    ))
    print("Practice_makes_perfect")
