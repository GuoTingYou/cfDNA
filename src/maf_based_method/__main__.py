import os
import sys
import time

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
