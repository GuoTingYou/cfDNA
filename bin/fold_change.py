import sys
import gzip
import numpy as np

def main():
    with gzip.open(sys.argv[1], 'rt') as f:
        f = (l.strip().split() for l in f)
        for chrom, start, end, *_, acc_p, acc_w in f:
            acc_p, acc_w = float(acc_p), float(acc_w)
            if acc_p == 0 or acc_w == 0:
                fc = 0
            else:
                fc = np.log(acc_w / acc_p) / np.log(2)
            print(chrom, start, end, acc_p, acc_w, fc, sep='\t')

if __name__ == "__main__":
    main()
