import sys

step = 20000
for chrom, length in map(lambda l: l.split(), sys.stdin):
    length = int(length)
    for start in range(0, length, step):
        end = start + step
        end = end if end < length else length
        print(chrom, start, end, sep='\t')
