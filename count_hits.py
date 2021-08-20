#!/home/max/bin/anaconda3/bin/python

import sys

entries = 0
nohits = 0

for line in sys.stdin:
    if line.startswith("SEQUENCE_ID="):
        entries +=1
    elif line.startswith("PRIMER_PAIR_NUM_RETURNED=0"):
        nohits += 1

print(f'{nohits}/{entries}')
